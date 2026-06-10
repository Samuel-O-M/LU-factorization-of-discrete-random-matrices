using Printf

const N = 4
const OUTPUT_FILE = "experiment/matrices_4x4.json"

function det_2x2(a11, a12, a21, a22)
    return a11 * a22 - a12 * a21
end

function det_3x3(m, i1, i2, i3, j1, j2, j3)
    return m[i1, j1] * det_2x2(m[i2, j2], m[i2, j3], m[i3, j2], m[i3, j3]) -
           m[i1, j2] * det_2x2(m[i2, j1], m[i2, j3], m[i3, j1], m[i3, j3]) +
           m[i1, j3] * det_2x2(m[i2, j1], m[i2, j2], m[i3, j1], m[i3, j2])
end

function det_4x4(m)
    return m[1,1] * det_3x3(m, 2,3,4, 2,3,4) -
           m[1,2] * det_3x3(m, 2,3,4, 1,3,4) +
           m[1,3] * det_3x3(m, 2,3,4, 1,2,4) -
           m[1,4] * det_3x3(m, 2,3,4, 1,2,3)
end

function is_strongly_nonsingular(mat::AbstractMatrix{Int})
    mat[1,1] != 1 && return false
    det_2x2(mat[1,1], mat[1,2], mat[2,1], mat[2,2]) == 0 && return false
    det_3x3(mat, 1,2,3, 1,2,3) == 0 && return false
    det_4x4(mat) == 0 && return false
    return true
end

function unpack!(mat::AbstractMatrix{Int}, p::UInt64)
    q = p
    for i in 1:N
        for j in 1:N
            mat[i, j] = q & 1
            q >>= 1
        end
    end
end

function write_json(matrices::Vector{Matrix{Int}})
    open(OUTPUT_FILE, "w") do io
        write(io, "[\n")
        for (idx, mat) in enumerate(matrices)
            write(io, "  [")
            for i in 1:N
                write(io, "[")
                for j in 1:N
                    write(io, string(mat[i, j]))
                    j < N && write(io, ",")
                end
                write(io, "]")
                i < N && write(io, ",")
            end
            write(io, "]")
            idx < length(matrices) && write(io, ",")
            write(io, "\n")
        end
        write(io, "]\n")
    end
end

function main()
    mat     = zeros(Int, N, N)
    results = Matrix{Int}[]

    total = UInt64(1) << (N * N - 1)
    for p_upper in UInt64(0):(total - 1)
        p = (p_upper << 1) | UInt64(1)
        unpack!(mat, p)
        is_strongly_nonsingular(mat) && push!(results, copy(mat))
    end

    @printf("Checked %d matrices, %d are strongly non‑singular  (%.4f %%)\n",
            total, length(results), 100 * length(results) / total)

    write_json(results)
    @printf("Written to %s\n", OUTPUT_FILE)
end

main()
