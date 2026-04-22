include("../utils/nauty.jl")
using .MyNauty
using StaticArrays
using Base.Threads
using Printf
using Dates
using SortingAlgorithms

const PARENTS_PER_BATCH = 2000
const MAGIC_FOOTER = 0xCAFEBABECAFEBABE
const NUM_SHARDS = 256

const LOCAL_SHARD_FLUSH_THRESH = 100_000
const GLOBAL_SHARD_COMPACT_THRESH = 20_000_000

const PackedMatrix = UInt64

struct Family
    rep::PackedMatrix
    count::UInt64
end

Base.write(io::IO, f::Family) = write(io, f.rep, f.count)

function Base.read(io::IO, ::Type{Family})
    r = read(io, PackedMatrix)
    c = read(io, UInt64)
    return Family(r, c)
end

Base.isless(a::Family, b::Family) = a.rep < b.rep

const ResultBuffer = Vector{Family}

function compact_results!(res::ResultBuffer)
    if isempty(res); return; end
    sort!(res, alg=RadixSort, by = x -> x.rep)

    write_idx = 1
    curr_rep = res[1].rep
    curr_count = res[1].count

    @inbounds for read_idx in 2:length(res)
        if res[read_idx].rep == curr_rep
            curr_count += res[read_idx].count
        else
            res[write_idx] = Family(curr_rep, curr_count)
            write_idx += 1
            curr_rep = res[read_idx].rep
            curr_count = res[read_idx].count
        end
    end

    @inbounds res[write_idx] = Family(curr_rep, curr_count)
    resize!(res, write_idx)
    return
end

@inline function unpack!(mat::AbstractMatrix, p::PackedMatrix, n::Int)
    mask = (UInt64(1) << n) - 1
    @inbounds for i in 0:n-1
        row_bits = (p >> (i * n)) & mask
        for j in 0:n-1
            mat[i+1, j+1] = Int32((row_bits >> j) & 1)
        end
    end
end

function det_bareiss_mut!(M::AbstractMatrix{Int64}, n::Int)
    if n == 1; return M[1, 1]; end
    if n == 2; return M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]; end
    if n == 3
        return M[1,1]*(M[2,2]*M[3,3] - M[2,3]*M[3,2]) -
               M[1,2]*(M[2,1]*M[3,3] - M[2,3]*M[3,1]) +
               M[1,3]*(M[2,1]*M[3,2] - M[2,2]*M[3,1])
    end

    sign = 1
    prev = Int64(1)
    @inbounds for k in 1:n-1
        if M[k, k] == 0
            pivot = 0
            for i in k+1:n
                if M[i, k] != 0
                    pivot = i
                    break
                end
            end
            if pivot == 0
                return Int64(0)
            end
            for j in k:n
                M[k, j], M[pivot, j] = M[pivot, j], M[k, j]
            end
            sign = -sign
        end

        p = M[k, k]
        for i in k+1:n
            for j in k+1:n
                M[i, j] = div(M[i, j] * p - M[i, k] * M[k, j], prev)
            end
        end
        prev = p
    end
    return sign * M[n, n]
end

function det_bareiss(mat::AbstractMatrix, n::Int, M::AbstractMatrix)
    @inbounds for j in 1:n, i in 1:n
        M[i, j] = mat[i, j]
    end
    return det_bareiss_mut!(M, n)
end

function faddeev_leverrier_adj!(mat::AbstractMatrix, n::Int, B::AbstractMatrix, M::AbstractMatrix, adjOut::AbstractMatrix)
    fill!(B, 0)
    @inbounds for i in 1:n; B[i, i] = 1; end
    @inbounds for k in 1:n
        for i in 1:n, j in 1:n
            s = Int64(0)
            for x in 1:n; s += Int64(mat[i, x]) * B[x, j]; end
            M[i, j] = s
        end

        tr = Int64(0)
        for i in 1:n; tr += M[i, i]; end
        ck = -div(tr, k)

        if k == n
            sv = ((n + 1) % 2 != 0) ? -1 : 1
            for i in 1:n, j in 1:n; adjOut[i, j] = sv * B[i, j]; end
        else
            for i in 1:n, j in 1:n
                v = M[i, j]
                if i == j; v += ck; end
                B[i, j] = v
            end
        end
    end
end


mutable struct ScratchPad
    mat::MMatrix{8, 8, Int32, 64}
    adjA::MMatrix{8, 8, Int64, 64}
    w::MVector{8, Int64}
    B::MMatrix{8, 8, Int64, 64}
    M::MMatrix{8, 8, Int64, 64}
    s1::MyNauty.Solver
    s2::MyNauty.Solver
    children_buffer::Vector{UInt64}

    function ScratchPad()
        new(@MMatrix(zeros(Int32, 8, 8)), @MMatrix(zeros(Int64, 8, 8)),
            @MVector(zeros(Int64, 8)), @MMatrix(zeros(Int64, 8, 8)),
            @MMatrix(zeros(Int64, 8, 8)),
            MyNauty.Solver(8), MyNauty.Solver(8),
            Vector{UInt64}(undef, 32768))
    end
end


function worker_compute!(parents::Vector{Family}, pn::Int, n::Int, local_res::ResultBuffer, pad::ScratchPad)
    combs = 1 << pn
    empty!(local_res)

    parent_counter = 0

    for fam in parents
        parent_counter += 1

        unpack!(pad.mat, fam.rep, pn)
        dA = det_bareiss(pad.mat, pn, pad.M)
        faddeev_leverrier_adj!(pad.mat, pn, pad.B, pad.M, pad.adjA)

        p_base = UInt64(0)
        @inbounds for i in 0:pn-1, j in 0:pn-1
            if pad.mat[i+1, j+1] != 0; p_base |= (UInt64(1) << (i * n + j)); end
        end

        v = UInt64(0)
        v_child_bits = UInt64(0)
        fill!(pad.w, 0)

        child_count = 0

        for kv in 0:(combs - 1)
            u = UInt64(0)
            u_child_bits = UInt64(0)
            S = Int64(0)

            v_child = p_base | v_child_bits

            for ku in 0:(combs - 1)
                child = v_child | u_child_bits

                c1_int = Int((S | -S) >>> 63)
                diff = dA - S
                c2_int = Int((diff | -diff) >>> 63)

                @inbounds pad.children_buffer[child_count + 1] = child
                child_count += c1_int

                @inbounds pad.children_buffer[child_count + 1] = child | (UInt64(1) << (pn * n + pn))
                child_count += c2_int

                if ku + 1 < combs
                    bit = trailing_zeros(ku + 1)
                    u_child_bits ⊻= (UInt64(1) << (pn * n + bit))
                    u ⊻= (UInt64(1) << bit)
                    S += ((u >> bit) & 1 != 0 ? 1 : -1) * pad.w[bit+1]
                end
            end

            if kv + 1 < combs
                bit = trailing_zeros(kv + 1)
                v ⊻= (UInt64(1) << bit)
                v_child_bits ⊻= (UInt64(1) << (bit * n + pn))

                sgn = ((v >> bit) & 1 != 0 ? 1 : -1)
                @inbounds for i in 1:pn
                    pad.w[i] += sgn * pad.adjA[i, bit+1]
                end
            end
        end

        @inbounds for i in 1:child_count
            canon = MyNauty.get_canonical!(pad.s1, pad.s2, pad.children_buffer[i], n)
            push!(local_res, Family(canon, fam.count))
        end

        if parent_counter % 100 == 0
            compact_results!(local_res)
        end
    end
end

function is_batch_valid(filename::String)
    if !isfile(filename); return false; end
    sz = stat(filename).size
    if sz < sizeof(UInt64); return false; end
    open(filename, "r") do io
        seek(io, sz - sizeof(UInt64))
        return read(io, UInt64) == MAGIC_FOOTER
    end
end

function log_results(n::Int)
    fname = "data/super/binary/families_$(n)x$(n).bin"
    if !is_batch_valid(fname); return; end

    sz = stat(fname).size
    unique_families = div(sz - sizeof(UInt64), sizeof(Family))
    total_matrices = UInt64(0)

    open(fname, "r") do io
        processed = 0
        buf_size = 10000
        buf = Vector{Family}(undef, buf_size)
        while processed < unique_families
            chunk = min(buf_size, unique_families - processed)
            GC.@preserve buf begin
                unsafe_read(io, pointer(buf), chunk * sizeof(Family))
            end
            @inbounds for i in 1:chunk
                total_matrices += buf[i].count
            end
            processed += chunk
        end
    end

    println(" [STATS] N=$n Families: $unique_families | Total Matrices: $total_matrices")
    open("results/super/summary.txt", "a") do log
        println(log, "N=$n Families=$unique_families Total=$total_matrices")
    end
end

function process_n(n::Int)
    mkpath("data/super/binary")
    pn = n - 1
    parent_file = "data/super/binary/families_$(pn)x$(pn).bin"
    final_target = "data/super/binary/families_$(n)x$(n).bin"

    if is_batch_valid(final_target)
        println("N=$n is already completed and verified.")
        log_results(n)
        return
    end

    if !isfile(parent_file)
        println(stderr, "Error: $parent_file not found")
        exit(1)
    end

    file_size = stat(parent_file).size
    total_families = div(file_size - sizeof(UInt64), sizeof(Family))
    total_batches = div(total_families + PARENTS_PER_BATCH - 1, PARENTS_PER_BATCH)

    println("=== Processing N=$n ===")
    println("Parents: $total_families | Estimated Batches: $total_batches")

    num_threads = Threads.nthreads()
    worker_slots = max(1, num_threads - 1)
    println("Threads: $num_threads | Workers: $worker_slots")

    global_shards = [ResultBuffer() for _ in 1:NUM_SHARDS]
    for k in 1:NUM_SHARDS; sizehint!(global_shards[k], GLOBAL_SHARD_COMPACT_THRESH + LOCAL_SHARD_FLUSH_THRESH); end
    shard_locks = [ReentrantLock() for _ in 1:NUM_SHARDS]
    job_channel = Channel{Vector{Family}}(worker_slots * 4)

    batches_processed = Atomic{Int}(0)
    start_time = now()

    worker_tasks =[]
    for _ in 1:worker_slots
        t = Threads.@spawn begin
            pad = ScratchPad()
            children = ResultBuffer()
            sizehint!(children, PARENTS_PER_BATCH * 1000)

            local_shards =[ResultBuffer() for _ in 1:NUM_SHARDS]
            for k in 1:NUM_SHARDS; sizehint!(local_shards[k], LOCAL_SHARD_FLUSH_THRESH); end

            for batch in job_channel
                worker_compute!(batch, pn, n, children, pad)
                compact_results!(children)

                shift_amt = max(0, n * n - 8)
                for child in children
                    idx = Int((child.rep >> shift_amt) & 0xFF) + 1

                    if length(local_shards[idx]) == LOCAL_SHARD_FLUSH_THRESH
                        lock(shard_locks[idx]) do
                            if length(global_shards[idx]) + length(local_shards[idx]) > GLOBAL_SHARD_COMPACT_THRESH
                                compact_results!(global_shards[idx])
                            end
                            append!(global_shards[idx], local_shards[idx])
                            if length(global_shards[idx]) > GLOBAL_SHARD_COMPACT_THRESH
                                compact_results!(global_shards[idx])
                            end
                        end
                        empty!(local_shards[idx])
                    end

                    push!(local_shards[idx], child)
                end

                bp = atomic_add!(batches_processed, 1) + 1
                if bp % 50 == 0 || bp == total_batches
                    elapsed = max((now() - start_time).value / 1000.0, 0.001)
                    speed = bp / elapsed
                    eta = (total_batches - bp) / speed
                    @printf("[Generation] %d/%d batches (%.1f%%) | ETA: %ds\n", bp, total_batches, (bp*100.0/total_batches), eta)
                    flush(stdout)
                end
            end

            for k in 1:NUM_SHARDS
                if !isempty(local_shards[k])
                    lock(shard_locks[k]) do
                        if length(global_shards[k]) + length(local_shards[k]) > GLOBAL_SHARD_COMPACT_THRESH
                            compact_results!(global_shards[k])
                        end
                        append!(global_shards[k], local_shards[k])
                        if length(global_shards[k]) > GLOBAL_SHARD_COMPACT_THRESH
                            compact_results!(global_shards[k])
                        end
                    end
                    empty!(local_shards[k])
                end
            end
        end
        push!(worker_tasks, t)
    end

    open(parent_file, "r") do in_stream
        bytes_total = file_size - sizeof(UInt64)
        while position(in_stream) < bytes_total
            bytes_left = bytes_total - position(in_stream)
            read_count = min(PARENTS_PER_BATCH, bytes_left ÷ sizeof(Family))

            parents = Vector{Family}(undef, read_count)
            GC.@preserve parents begin
                unsafe_read(in_stream, pointer(parents), read_count * sizeof(Family))
            end
            put!(job_channel, parents)
        end
    end

    close(job_channel)
    for t in worker_tasks; wait(t); end

    println("\n[Aggregator] Generation complete. Parallel sorting 256 memory shards (max 16 concurrent)...")

    sem = Base.Semaphore(16)
    sort_tasks =[]
    for k in 1:NUM_SHARDS
        t = Threads.@spawn begin
            Base.acquire(sem)
            try
                compact_results!(global_shards[k])
            finally
                Base.release(sem)
            end
        end
        push!(sort_tasks, t)
    end
    for t in sort_tasks; wait(t); end

    println("[Aggregator] Sort complete. Assembling final binary file directly from RAM...")
    tmp_target = final_target * ".tmp"
    open(tmp_target, "w") do out
        for k in 1:NUM_SHARDS
            shard = global_shards[k]
            if !isempty(shard)
                GC.@preserve shard begin
                    unsafe_write(out, pointer(shard), sizeof(Family) * length(shard))
                end
            end
            empty!(shard)
            sizehint!(shard, 0)
            GC.gc(false)
        end
        write(out, MAGIC_FOOTER)
    end

    mv(tmp_target, final_target; force=true)
    println("Process N=$n FINISHED. Result in $final_target")
    log_results(n)
end

function main()
    mkpath("data/super/binary")
    mkpath("results/super")

    f1 = "data/super/binary/families_1x1.bin"
    if !is_batch_valid(f1)
        f = Family(UInt64(1), UInt64(1))
        tmp_name = f1 * ".tmp"
        open(tmp_name, "w") do out
            write(out, f)
            write(out, MAGIC_FOOTER)
        end
        mv(tmp_name, f1; force=true)
        println("Initializing N=1")
        open("results/super/summary.txt", "w") do log
            println(log, "N=1 Families=1 Total=1")
        end
    else
        open("results/super/summary.txt", "w") do log
            println(log, "N=1 Families=1 Total=1")
        end
        log_results(1)
    end

    for n in 2:8
        process_n(n)
    end
end

main()
