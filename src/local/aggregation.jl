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

struct Signature
    r::SVector{8, UInt8}
    c::SVector{8, UInt8}
end
Base.isless(a::Signature, b::Signature) = a.r != b.r ? a.r < b.r : a.c < b.c

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

function aggregate_n8(parents::Vector{Family})
    filename = "data/local/binary/counts_8x8.bin"
    pn, combs = 7, 1 << 7
    num_slots = Threads.maxthreadid()
    local_maps = [Dict{Signature, UInt64}() for _ in 1:num_slots]
    pads = [ScratchPad() for _ in 1:num_slots]
    
    @threads :static for p_idx in 1:length(parents)
        tid = Threads.threadid()
        pad, lmap, fam = pads[tid], local_maps[tid], parents[p_idx]
        unpack!(pad.mat, fam.rep, pn)
        dA = det_fast(pad.mat, pn, pad.det_stack)
        faddeev_leverrier_adj!(pad.mat, pn, pad)
        fill!(pad.base_r, 0); fill!(pad.base_c, 0)
        for i in 1:pn, j in 1:pn
            if pad.mat[i, j] != 0; pad.base_r[i] += 1; pad.base_c[j] += 1; end
        end
        v, w, current_c7 = UInt64(0), pad.w, 0
        fill!(w, 0); pad.r .= pad.base_r
        for kv in 0:(combs - 1)
            u, S, current_r7 = UInt64(0), Int64(0), 0
            pad.c .= pad.base_c
            for ku in 0:(combs - 1)
                if S != 0
                    pad.r[8], pad.c[8] = current_r7, current_c7
                    sig = Signature(sort(SVector{8, UInt8}(pad.r)), sort(SVector{8, UInt8}(pad.c)))
                    lmap[sig] = get(lmap, sig, 0) + fam.count
                end
                if (dA - S) != 0
                    pad.r[8], pad.c[8] = current_r7 + 1, current_c7 + 1
                    sig = Signature(sort(SVector{8, UInt8}(pad.r)), sort(SVector{8, UInt8}(pad.c)))
                    lmap[sig] = get(lmap, sig, 0) + fam.count
                end
                if ku + 1 < combs
                    bit = trailing_zeros(ku + 1); u ⊻= (UInt64(1) << bit)
                    sgn = ((u >> bit) & 1 != 0 ? 1 : -1)
                    S += sgn * w[bit+1]; pad.c[bit+1] += sgn; current_r7 += sgn
                end
            end
            if kv + 1 < combs
                bit = trailing_zeros(kv + 1); v ⊻= (UInt64(1) << bit)
                sgn = ((v >> bit) & 1 != 0 ? 1 : -1)
                for i in 1:pn; w[i] += sgn * pad.adjA[i, bit+1]; end
                pad.r[bit+1] += sgn; current_c7 += sgn
            end
        end
    end
    gm = Dict{Signature, UInt64}()
    for lm in local_maps, (sig, count) in lm; gm[sig] = get(gm, sig, 0) + count; end
    sorted_sigs = sort(collect(keys(gm)))
    total = UInt64(0)
    open(filename, "w") do io
        for sig in sorted_sigs; c = gm[sig]; total += c; write(io, sig.r); write(io, sig.c); write(io, c); end
    end
    num_signatures = length(gm)
    return (total, num_signatures)
end

function main()
    mkpath("data/local/binary")
    mkpath("results/local")
    
    parents = load_families("data/local/binary/families_7x7.bin")
    if isempty(parents)
        println(stderr, "Error: data/local/binary/families_7x7.bin not found. Run families.jl first.")
        exit(1)
    end
    
    total_8, num_signatures = aggregate_n8(parents)
    
    open("results/local/summary.txt", "a") do report
        println(report, "N=8 Signatures=$num_signatures Total=$total_8")
    end
    @printf("N=8 Signatures=%d Total=%d\n", num_signatures, total_8)
end

main()
