using StaticArrays
using Base.Threads
using Printf
using Dates

const PARENTS_PER_BATCH = 5000
const MAGIC_FOOTER = 0xCAFEBABECAFEBABE
const NUM_SHARDS = 256

const PackedMatrix = UInt64

struct Family
    rep::PackedMatrix
    count::UInt64
end

struct Signature
    r::SVector{9, UInt8}
    c::SVector{9, UInt8}
end

Base.isless(a::Signature, b::Signature) = a.r != b.r ? a.r < b.r : a.c < b.c
Base.hash(s::Signature, h::UInt) = hash(s.c, hash(s.r, h))
Base.:(==)(a::Signature, b::Signature) = a.r == b.r && a.c == b.c

@inline function unpack!(mat::AbstractMatrix{Int64}, p::PackedMatrix, n::Int)
    mask = (UInt64(1) << n) - 1
    @inbounds for i in 0:n-1
        row_bits = (p >> (i * n)) & mask
        for j in 0:n-1
            mat[i+1, j+1] = Int64((row_bits >> j) & 1)
        end
    end
end

function det_bareiss_mut!(M::AbstractMatrix{Int64}, n::Int)::Int64
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
            if pivot == 0; return Int64(0); end
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

function det_bareiss(mat::AbstractMatrix{Int64}, n::Int, M::AbstractMatrix{Int64})::Int64
    @inbounds for j in 1:n, i in 1:n
        M[i, j] = mat[i, j]
    end
    return det_bareiss_mut!(M, n)
end

function faddeev_leverrier_adj!(mat::AbstractMatrix{Int64}, n::Int, B::AbstractMatrix{Int64}, M::AbstractMatrix{Int64}, adjOut::AbstractMatrix{Int64})
    fill!(B, 0)
    @inbounds for i in 1:n; B[i, i] = 1; end
    @inbounds for k in 1:n
        for i in 1:n, j in 1:n
            s = Int64(0)
            for x in 1:n; s += mat[i, x] * B[x, j]; end
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
    mat::MMatrix{8, 8, Int64, 64}
    adjA::MMatrix{8, 8, Int64, 64}
    w::MVector{8, Int64}
    B::MMatrix{8, 8, Int64, 64}
    M::MMatrix{8, 8, Int64, 64}
    r::MVector{9, UInt8}
    c::MVector{9, UInt8}
    base_r::MVector{8, UInt8}
    base_c::MVector{8, UInt8}

    function ScratchPad()
        new(@MMatrix(zeros(Int64, 8, 8)), @MMatrix(zeros(Int64, 8, 8)),
            @MVector(zeros(Int64, 8)), @MMatrix(zeros(Int64, 8, 8)),
            @MMatrix(zeros(Int64, 8, 8)),
            @MVector(zeros(UInt8, 9)), @MVector(zeros(UInt8, 9)),
            @MVector(zeros(UInt8, 8)), @MVector(zeros(UInt8, 8)))
    end
end

function worker_aggregate!(parents::Vector{Family}, pn::Int, local_map::Dict{Signature, UInt128}, pad::ScratchPad)
    combs = 1 << pn
    empty!(local_map)

    for fam in parents
        unpack!(pad.mat, fam.rep, pn)
        dA = det_bareiss(pad.mat, pn, pad.M)
        faddeev_leverrier_adj!(pad.mat, pn, pad.B, pad.M, pad.adjA)

        fill!(pad.base_r, 0); fill!(pad.base_c, 0)
        @inbounds for i in 1:pn, j in 1:pn
            if pad.mat[i, j] != 0
                pad.base_r[i] += 1
                pad.base_c[j] += 1
            end
        end

        v = UInt64(0)
        current_c9 = UInt8(0)
        fill!(pad.w, 0)
        pad.r[1:8] .= pad.base_r

        count_u128 = UInt128(fam.count)

        for kv in 0:(combs - 1)
            u = UInt64(0); S = Int64(0); current_r9 = UInt8(0)
            pad.c[1:8] .= pad.base_c

            for ku in 0:(combs - 1)
                if S != 0
                    pad.r[9] = current_r9
                    pad.c[9] = current_c9
                    sig = Signature(sort(SVector{9, UInt8}(pad.r)), sort(SVector{9, UInt8}(pad.c)))
                    local_map[sig] = get(local_map, sig, zero(UInt128)) + count_u128
                end
                if (dA - S) != 0
                    pad.r[9] = current_r9 + 1
                    pad.c[9] = current_c9 + 1
                    sig = Signature(sort(SVector{9, UInt8}(pad.r)), sort(SVector{9, UInt8}(pad.c)))
                    local_map[sig] = get(local_map, sig, zero(UInt128)) + count_u128
                end

                if ku + 1 < combs
                    bit = trailing_zeros(ku + 1)
                    u ⊻= (UInt64(1) << bit)
                    sgn = ((u >> bit) & 1 != 0 ? 1 : -1)
                    S += sgn * pad.w[bit+1]
                    pad.c[bit+1] = UInt8(Int(pad.c[bit+1]) + sgn)
                    current_r9 = UInt8(Int(current_r9) + sgn)
                end
            end

            if kv + 1 < combs
                bit = trailing_zeros(kv + 1)
                v ⊻= (UInt64(1) << bit)
                sgn = ((v >> bit) & 1 != 0 ? 1 : -1)
                @inbounds for i in 1:pn
                    pad.w[i] += sgn * pad.adjA[i, bit+1]
                end
                pad.r[bit+1] = UInt8(Int(pad.r[bit+1]) + sgn)
                current_c9 = UInt8(Int(current_c9) + sgn)
            end
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

function main()
    mkpath("data/super/binary")
    mkpath("results/super")
    parent_file = "data/super/binary/families_8x8.bin"
    output_file = "data/super/binary/counts_9x9.bin"

    if !is_batch_valid(parent_file)
        println(stderr, "Error: $parent_file not found or incomplete.")
        exit(1)
    end

    file_size = stat(parent_file).size
    total_families = div(file_size - sizeof(UInt64), sizeof(Family))
    total_batches = div(total_families + PARENTS_PER_BATCH - 1, PARENTS_PER_BATCH)

    println("=== Commencing 9x9 Signature Aggregation ===")
    println("Parents: $total_families | Batches: $total_batches")

    num_threads = Threads.nthreads()
    worker_slots = max(1, num_threads - 1)
    println("Threads: $num_threads | Workers: $worker_slots")

    global_maps = [Dict{Signature, UInt128}() for _ in 1:NUM_SHARDS]
    shard_locks = [ReentrantLock() for _ in 1:NUM_SHARDS]
    job_channel = Channel{Vector{Family}}(worker_slots * 4)

    batches_processed = Atomic{Int}(0)
    start_time = now()

    worker_tasks = []
    for _ in 1:worker_slots
        t = Threads.@spawn begin
            pad = ScratchPad()
            local_map = Dict{Signature, UInt128}()
            sizehint!(local_map, 200_000)

            for batch in job_channel
                worker_aggregate!(batch, 8, local_map, pad)

                for (sig, count) in local_map
                    shard_idx = (hash(sig) % NUM_SHARDS) + 1
                    lock(shard_locks[shard_idx]) do
                        global_maps[shard_idx][sig] = get(global_maps[shard_idx], sig, zero(UInt128)) + count
                    end
                end

                bp = atomic_add!(batches_processed, 1) + 1
                if (bp % (bp <= 1000 ? 50 : 1000) == 0 || bp == total_batches)
                    elapsed = max((now() - start_time).value / 1000.0, 0.001)
                    speed = bp / elapsed
                    eta = (total_batches - bp) / speed
                    @printf("[Aggregation] %d/%d batches (%.1f%%) | ETA: %ds\n", bp, total_batches, (bp*100.0/total_batches), eta)
                    flush(stdout)
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

    println("\n[IO Phase] Aggregation complete. Writing unified UInt128 records to disk...")
    total_matrices = UInt128(0)

    tmp_output = output_file * ".tmp"
    open(tmp_output, "w") do out
        for k in 1:NUM_SHARDS
            shard = global_maps[k]
            for sig in sort!(collect(keys(shard)))
                c = shard[sig]
                total_matrices += c
                write(out, sig.r)
                write(out, sig.c)
                write(out, c)
            end
            empty!(shard)
            GC.gc(false)
        end
    end

    num_signatures = sum(length(shard) for shard in global_maps)
    mv(tmp_output, output_file; force=true)
    println("Process FINISHED. Result in $output_file")
    @printf("N=9 Signatures=%d Total=%s\n", num_signatures, string(total_matrices))
    open("results/super/summary.txt", "a") do log
        println(log, "N=9 Signatures=$num_signatures Total=$total_matrices")
    end
end

main()
