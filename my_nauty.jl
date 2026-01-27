module MyNauty

using StaticArrays

const MAX_V = 16

struct Partition
    lab::SVector{MAX_V, UInt8}
    ptn::SVector{MAX_V, UInt8}
end

mutable struct Solver
    N::Int
    V::Int
    adj::MVector{MAX_V, UInt16} 
    best_key::UInt64
    has_best::Bool
    best_part::Partition
    # orbits[depth, vertex]
    orbits::MMatrix{17, MAX_V, UInt8, 272}
    fixed_path::MVector{MAX_V, UInt8}
    
    # Buffers for refinement to avoid allocation
    vertex_weights::MVector{MAX_V, UInt64}
    sigs::MVector{MAX_V, UInt64}
    gamma::MVector{MAX_V, UInt8}

    function Solver(n::Int)
        v_size = 2 * n
        initial_lab = @SVector zeros(UInt8, MAX_V)
        initial_ptn = @SVector zeros(UInt8, MAX_V)
        new(n, v_size, @MVector(zeros(UInt16, MAX_V)), typemax(UInt64), false, 
            Partition(initial_lab, initial_ptn), @MMatrix(zeros(UInt8, 17, MAX_V)), 
            @MVector(zeros(UInt8, MAX_V)), @MVector(zeros(UInt64, MAX_V)), 
            @MVector(zeros(UInt64, MAX_V)), @MVector(zeros(UInt8, MAX_V)))
    end
end

function reset!(s::Solver, matrix_bits::UInt64, n::Int)
    s.N = n
    s.V = 2n
    s.has_best = false
    s.best_key = typemax(UInt64)
    fill!(s.adj, 0)
    @inbounds for r in 0:n-1
        for c in 0:n-1
            if ((matrix_bits >> (r * n + c)) & 1) != 0
                s.adj[r + 1] |= (UInt16(1) << (n + c))
                s.adj[(n + c) + 1] |= (UInt16(1) << r)
            end
        end
    end
    @inbounds for i in 0:MAX_V-1; s.orbits[1, i + 1] = UInt8(i); end
end

@inline function find_orbit!(s::Solver, k::Int, v::Integer)::UInt8
    d = k + 1
    root = UInt8(v)
    @inbounds while root != s.orbits[d, root + 1]; root = s.orbits[d, root + 1]; end
    curr = UInt8(v)
    @inbounds while curr != root
        nxt = s.orbits[d, curr + 1]
        s.orbits[d, curr + 1] = root
        curr = nxt
    end
    return root
end

@inline function join_orbits!(s::Solver, k::Int, u::Integer, v::Integer)
    ru, rv = find_orbit!(s, k, u), find_orbit!(s, k, v)
    if ru != rv
        d = k + 1
        @inbounds ru < rv ? (s.orbits[d, rv + 1] = ru) : (s.orbits[d, ru + 1] = rv)
    end
end

function update_automorphism!(s::Solver, leaf::Partition, depth::Int)
    @inbounds for i in 1:s.V; s.gamma[s.best_part.lab[i] + 1] = leaf.lab[i]; end
    @inbounds for k in 0:depth-1
        fixed_v = s.fixed_path[k + 1]
        if s.gamma[fixed_v + 1] != fixed_v; break; end
        for v in 0:s.V-1; join_orbits!(s, k + 1, v, s.gamma[v + 1]); end
    end
end

function refine!(s::Solver, p::Partition)
    lab, ptn = MVector(p.lab), MVector(p.ptn)
    changed = true
    @inbounds while changed
        changed = false
        cell_idx = 0
        for i in 1:s.V
            s.vertex_weights[lab[i] + 1] = UInt64(1) << (4 * cell_idx)
            if ptn[i] == 0; cell_idx += 1; end
        end
        ptr = 1
        while ptr <= s.V
            nxt = ptr
            while ptn[nxt] != 0; nxt += 1; end
            if ptr != nxt
                for i in ptr:nxt
                    u, sum_w, neighbors = lab[i], UInt64(0), s.adj[lab[i]+1]
                    while neighbors != 0
                        v_idx = trailing_zeros(neighbors)
                        sum_w += s.vertex_weights[v_idx + 1]
                        neighbors &= neighbors - UInt16(1)
                    end
                    s.sigs[i] = sum_w
                end
                for i in (ptr + 1):nxt # Insertion sort
                    v_l, v_s, j = lab[i], s.sigs[i], i - 1
                    while j >= ptr && s.sigs[j] > v_s
                        lab[j+1], s.sigs[j+1] = lab[j], s.sigs[j]
                        j -= 1
                    end
                    lab[j+1], s.sigs[j+1] = v_l, v_s
                end
                for i in ptr:(nxt - 1)
                    if s.sigs[i] != s.sigs[i + 1]; ptn[i] = 0; changed = true; end
                end
            end
            ptr = nxt + 1
        end
    end
    is_discrete = true
    @inbounds for i in 1:s.V; if ptn[i] != 0; is_discrete = false; break; end; end
    return is_discrete, Partition(SVector(lab), SVector(ptn))
end

function search!(s::Solver, p::Partition, depth::Int)
    is_discrete, p_ref = refine!(s, p)
    if is_discrete
        key = UInt64(0)
        @inbounds for r in 0:s.N-1
            u, neighbors = p_ref.lab[r + 1], s.adj[p_ref.lab[r + 1] + 1]
            for c in 0:s.N-1
                if ((neighbors >> p_ref.lab[(s.N + c) + 1]) & 1) != 0
                    key |= (UInt64(1) << (r * s.N + c))
                end
            end
        end
        if !s.has_best || key < s.best_key
            s.best_key, s.best_part, s.has_best = key, p_ref, true
        elseif key == s.best_key
            update_automorphism!(s, p_ref, depth)
        end
        return
    end
    t_start, t_end, min_sz, ptr = -1, -1, 99, 1
    @inbounds while ptr <= s.V
        nxt = ptr
        while p_ref.ptn[nxt] != 0; nxt += 1; end
        sz = nxt - ptr + 1
        if sz > 1 && sz < min_sz; min_sz, t_start, t_end = sz, ptr, nxt; end
        if sz == 2; break; end
        ptr = nxt + 1
    end
    if t_start == -1; return; end
    @inbounds s.orbits[depth + 2, :] .= s.orbits[depth + 1, :]
    @inbounds for i in t_start:t_end
        v = p_ref.lab[i]
        if find_orbit!(s, depth + 1, v) != v; continue; end
        c_lab, c_ptn = MVector(p_ref.lab), MVector(p_ref.ptn)
        c_lab[t_start], c_lab[i] = c_lab[i], c_lab[t_start]
        c_ptn[t_start] = 0
        s.fixed_path[depth + 1] = v
        search!(s, Partition(SVector(c_lab), SVector(c_ptn)), depth + 1)
    end
end

function get_canonical!(s1::Solver, s2::Solver, bits::UInt64, n::Int)
    reset!(s1, bits, n)
    l1, p1 = MVector{MAX_V, UInt8}(undef), MVector{MAX_V, UInt8}(undef)
    for i in 0:2n-1; l1[i+1], p1[i+1] = i, 1; end
    p1[n], p1[2n] = 0, 0
    search!(s1, Partition(SVector(l1), SVector(p1)), 0)
    
    reset!(s2, bits, n)
    l2, p2 = MVector{MAX_V, UInt8}(undef), MVector{MAX_V, UInt8}(undef)
    for i in 0:n-1; l2[i+1], l2[n+i+1] = n+i, i; p2[i+1], p2[n+i+1] = 1, 1; end
    p2[n], p2[2n] = 0, 0
    search!(s2, Partition(SVector(l2), SVector(p2)), 0)
    return min(s1.best_key, s2.best_key)
end

end