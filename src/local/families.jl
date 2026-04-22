include("../utils/nauty.jl")
using .MyNauty
using StaticArrays
using Base.Threads
using Printf

struct Family
    rep::UInt64
    count::UInt64
end
Base.isless(a::Family, b::Family) = a.rep < b.rep

mutable struct ScratchPad
    mat::MMatrix{8, 8, Int32, 64}
    adjA::MMatrix{8, 8, Int64, 64}
    w::MVector{8, Int64} 
    B::MMatrix{8, 8, Int64, 64}
    M::MMatrix{8, 8, Int64, 64}
    r::MVector{8, Int32}
    c::MVector{8, Int32}
    base_r::MVector{8, Int32}
    base_c::MVector{8, Int32}
    s1::MyNauty.Solver
    s2::MyNauty.Solver
    det_stack::Vector{MMatrix{8, 8, Int32, 64}}

    function ScratchPad()
        new(@MMatrix(zeros(Int32, 8, 8)), @MMatrix(zeros(Int64, 8, 8)),
            @MVector(zeros(Int64, 8)), @MMatrix(zeros(Int64, 8, 8)),
            @MMatrix(zeros(Int64, 8, 8)), @MVector(zeros(Int32, 8)),
            @MVector(zeros(Int32, 8)), @MVector(zeros(Int32, 8)),
            @MVector(zeros(Int32, 8)), MyNauty.Solver(8), MyNauty.Solver(8),
            [MMatrix{8, 8, Int32, 64}(undef) for _ in 1:8])
    end
end

function save_families(filename::String, families::Vector{Family})
    open(filename, "w") do io; write(io, families); end
end

function load_families(filename::String)::Vector{Family}
    if !isfile(filename); return Family[]; end
    fs = stat(filename).size
    fams = Vector{Family}(undef, div(fs, sizeof(Family)))
    open(filename, "r") do io; read!(io, fams); end
    return fams
end

@inline function unpack!(mat::AbstractMatrix, p::UInt64, n::Int)
    @inbounds for i in 0:n-1, j in 0:n-1
        mat[i+1, j+1] = Int32((p >> (i * n + j)) & 1)
    end
end

function det_fast(mat::AbstractMatrix, n::Int, stack::Vector{MMatrix{8, 8, Int32, 64}}, depth::Int=1)::Int64
    if n == 1; return @inbounds Int64(mat[1, 1]); end
    if n == 2; @inbounds return Int64(mat[1, 1]) * mat[2, 2] - Int64(mat[1, 2]) * mat[2, 1]; end
    d, sign, sub = Int64(0), 1, stack[depth]
    @inbounds for c in 1:n
        if mat[1, c] != 0
            for i in 2:n
                sc = 1
                for j in 1:n
                    if j == c; continue; end
                    sub[i-1, sc] = mat[i, j]
                    sc += 1
                end
            end
            term = Int64(mat[1, c]) * det_fast(sub, n - 1, stack, depth + 1)
            d += (sign == 1 ? term : -term)
        end
        sign = -sign
    end
    return d
end

function faddeev_leverrier_adj!(mat::AbstractMatrix, n::Int, pad::ScratchPad)
    B, M, adjOut = pad.B, pad.M, pad.adjA
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

function generate_level(parents::Vector{Family}, n::Int)
    mkpath("data/local/binary")
    filename = "data/local/binary/families_$(n)x$(n).bin"
    if isfile(filename); return load_families(filename); end
    
    pn, combs = n - 1, 1 << (n - 1)
    num_slots = Threads.maxthreadid()
    thread_results = [Vector{Family}() for _ in 1:num_slots]
    pads = [ScratchPad() for _ in 1:num_slots]
    
    println("Generating N=$n using $(Threads.nthreads()) threads...")
    @threads :static for p_idx in 1:length(parents)
        tid = Threads.threadid()
        pad, res, fam = pads[tid], thread_results[tid], parents[p_idx]
        unpack!(pad.mat, fam.rep, pn)
        dA = det_fast(pad.mat, pn, pad.det_stack)
        faddeev_leverrier_adj!(pad.mat, pn, pad)
        
        p_base = UInt64(0)
        for i in 0:pn-1, j in 0:pn-1
            if pad.mat[i+1, j+1] != 0; p_base |= (UInt64(1) << (i * n + j)); end
        end
        
        v, w = UInt64(0), pad.w
        fill!(w, 0)
        for kv in 0:(combs - 1)
            u, S = UInt64(0), Int64(0)
            for ku in 0:(combs - 1)
                if S != 0 || (dA - S) != 0
                    child = p_base
                    tv, tu = v, u
                    while tv != 0; b = trailing_zeros(tv); child |= (UInt64(1) << (b * n + pn)); tv &= tv - 1; end
                    while tu != 0; b = trailing_zeros(tu); child |= (UInt64(1) << (pn * n + b)); tu &= tu - 1; end
                    if S != 0
                        push!(res, Family(MyNauty.get_canonical!(pad.s1, pad.s2, child, n), fam.count))
                    end
                    if (dA - S) != 0
                        push!(res, Family(MyNauty.get_canonical!(pad.s1, pad.s2, child | (UInt64(1) << (pn * n + pn)), n), fam.count))
                    end
                end
                if ku + 1 < combs
                    bit = trailing_zeros(ku + 1); u ⊻= (UInt64(1) << bit)
                    S += ((u >> bit) & 1 != 0 ? 1 : -1) * w[bit+1]
                end
            end
            if kv + 1 < combs
                bit = trailing_zeros(kv + 1); v ⊻= (UInt64(1) << bit)
                sgn = ((v >> bit) & 1 != 0 ? 1 : -1)
                for i in 1:pn; w[i] += sgn * pad.adjA[i, bit+1]; end
            end
        end
    end
    
    println("Merging thread results...")
    all_children = Vector{Family}()
    sizehint!(all_children, sum(length, thread_results; init=0))
    for t_res in thread_results; append!(all_children, t_res); end
    sort!(all_children, by = x -> x.rep)
    
    unique_families = Vector{Family}()
    if !isempty(all_children)
        cr, cc = all_children[1].rep, all_children[1].count
        for i in 2:length(all_children)
            if all_children[i].rep == cr; cc += all_children[i].count
            else; push!(unique_families, Family(cr, cc)); cr, cc = all_children[i].rep, all_children[i].count; end
        end
        push!(unique_families, Family(cr, cc))
    end
    save_families(filename, unique_families)
    return unique_families
end

function main()
    mkpath("data/local/binary")
    mkpath("results/local")
    f1 = "data/local/binary/families_1x1.bin"
    current = [Family(UInt64(1), UInt64(1))]
    save_families(f1, current)
    total = sum(f.count for f in current; init=UInt64(0))
    open("results/local/summary.txt", "w") do report
        println(report, "N=1 Families=$(length(current)) Total=$total")
    end
    @printf("N=1 Families=%d Total=%d\n", length(current), total)
    for n in 2:7
        current = generate_level(current, n)
        total = sum(f.count for f in current; init=UInt64(0))
        open("results/local/summary.txt", "a") do report
            println(report, "N=$n Families=$(length(current)) Total=$total")
        end
        @printf("N=%d Families=%d Total=%d\n", n, length(current), total)
    end
end

main()
