using Base.Threads
using Printf
using Random

const TARGET_N      = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 30
const X_MULTIPLIER  = 200.0
const BASE_N        = 4
const DB_FILE       = "experiment/matrices_4x4.json"
const OUTPUT_FILE   = "experiment/ber_n$(TARGET_N).json"
const GLOBAL_SEED   = 42

const MILESTONES = Int[6, 12, 20, 30, 45, 70]

function expansion_sizes(target::Int)
    sizes = Int[BASE_N]
    for m in MILESTONES
        m > sizes[end] && m <= target && push!(sizes, m)
    end
    sizes[end] != target && push!(sizes, target)
    return sizes[2:end]
end

const BIG_ZERO = BigInt(0)
const BIG_ONE  = BigInt(1)

function load_4x4_matrices(filepath::String)
    text  = read(filepath, String)
    clean = replace(text, r"\s+" => "")
    row_matches = eachmatch(r"\[(0|1),(0|1),(0|1),(0|1)\]", clean)
    row_vals    = [parse.(Int, split(m.match[2:end-1], ',')) for m in row_matches]
    n           = length(row_vals) ÷ 4
    matrices    = Matrix{Int}[]
    for i in 1:n
        mat = Matrix{Int}(undef, 4, 4)
        for r in 1:4
            vals = row_vals[(i-1)*4 + r]
            for c in 1:4
                mat[r, c] = vals[c]
            end
        end
        push!(matrices, mat)
    end
    return matrices
end

function first_singular_minor(A::Matrix{BigInt}, sz::Int)
    sz == 0 && return 0
    temp       = A[1:sz, 1:sz]
    prev_pivot = BigInt(1)
    for k in 1:sz
        current_pivot = temp[k, k]
        iszero(current_pivot) && return k
        if k < sz
            for i in k+1:sz
                for j in k+1:sz
                    temp[i, j] = (temp[k, k] * temp[i, j] -
                                  temp[i, k] * temp[k, j]) ÷ prev_pivot
                end
            end
            prev_pivot = current_pivot
        end
    end
    return 0
end

function run_single_p(matrices::Vector{Matrix{Int}}, ones_counts::Vector{Int},
                      p::BigFloat, check_sizes::Vector{Int}, target_n::Int,
                      x_mult::Float64,
                      rngs::Vector{Xoshiro}, buffers::Vector{Matrix{BigInt}})
    nslots     = Threads.maxthreadid()
    n_matrices = length(matrices)

    p_float = Float64(p)

    main_rng = rngs[1]

    workloads = Tuple{Int, Int}[]
    total_ext = 0
    for i in 1:n_matrices
        expected = x_mult
        n_ext    = floor(Int, expected)
        rand(main_rng) < (expected - n_ext) && (n_ext += 1)
        if n_ext > 0
            push!(workloads, (i, n_ext))
            total_ext += n_ext
        end
    end
    shuffle!(main_rng, workloads)

    thread_success   = zeros(Int, nslots)
    thread_prob_N    = zeros(BigFloat, nslots)
    thread_fail_by_k = [zeros(Int, target_n) for _ in 1:nslots]

    Threads.@threads for job in workloads
        tid       = Threads.threadid()
        rng       = rngs[tid]
        A         = buffers[tid]
        mat_idx   = job[1]
        n_ext     = job[2]

        base = matrices[mat_idx]

        o = ones_counts[mat_idx]
        prob_A = (p^o) * ((BigFloat(1.0) - p)^(16 - o))
        weight = prob_A / BigFloat(x_mult)

        local_success = 0

        for _ in 1:n_ext
            for r in 1:BASE_N, c in 1:BASE_N
                A[r, c] = base[r, c] == 1 ? BIG_ONE : BIG_ZERO
            end

            prev_size = BASE_N
            k         = 0
            for cur_size in check_sizes
                for r in prev_size+1:cur_size
                    for c in 1:cur_size
                        A[r, c] = rand(rng) < p_float ? BIG_ONE : BIG_ZERO
                    end
                end
                for r in 1:prev_size
                    for c in prev_size+1:cur_size
                        A[r, c] = rand(rng) < p_float ? BIG_ONE : BIG_ZERO
                    end
                end

                k = first_singular_minor(A, cur_size)
                k != 0 && break
                prev_size = cur_size
            end

            if k == 0
                local_success += 1
            else
                thread_fail_by_k[tid][k] += 1
            end
        end

        thread_success[tid] += local_success
        thread_prob_N[tid] += local_success * weight
    end

    total_success = sum(thread_success)
    total_failed  = total_ext - total_success

    fail_by_k = zeros(Int, target_n)
    for t in 1:nslots
        fail_by_k .+= thread_fail_by_k[t]
    end

    p_4 = BigFloat(0)
    for i in 1:n_matrices
        o = ones_counts[i]
        p_4 += (p^o) * ((BigFloat(1.0) - p)^(16 - o))
    end

    p_n = sum(thread_prob_N)

    p_cond = p_4 > 0.0 ? p_n / p_4 : BigFloat(0)

    return p_n, p_4, p_cond, total_success, total_ext, fail_by_k
end

function write_json(results::Vector{Dict{Symbol, Any}})
    open(OUTPUT_FILE, "w") do io
        write(io, "{\n")
        @printf(io, "  \"parameters\": {\n")
        @printf(io, "    \"target_n\": %d,\n", TARGET_N)
        @printf(io, "    \"x_multiplier\": %.2f,\n", X_MULTIPLIER)
        @printf(io, "    \"base_n\": %d,\n", BASE_N)
        @printf(io, "    \"seed\": %d\n", GLOBAL_SEED)
        write(io, "  },\n")
        write(io, "  \"results\": [\n")

        for (idx, res) in enumerate(results)
            write(io, "    {\n")
            @printf(io, "      \"p\": %.2f,\n",          Float64(res[:p]))
            @printf(io, "      \"probability\": %s,\n",  string(res[:probability]))
            @printf(io, "      \"p4_sns\": %s,\n",       string(res[:p4_sns]))
            @printf(io, "      \"conditional\": %s,\n",  string(res[:conditional]))
            @printf(io, "      \"extensions\": %d,\n",   res[:extensions])
            @printf(io, "      \"successes\": %d\n",     res[:successes])
            write(io, "    }")
            idx < length(results) && write(io, ",")
            write(io, "\n")
        end

        write(io, "  ]\n")
        write(io, "}\n")
    end
    println("[ber] Results saved to $OUTPUT_FILE")
end

function main()
    setprecision(BigFloat, 512)

    P_VALUES = [BigFloat(i) / BigFloat(100) for i in 1:100]

    matrices = load_4x4_matrices(DB_FILE)
    n_matrices = length(matrices)
    @printf("[ber] Loaded %d strongly non‑singular 4×4 matrices\n", n_matrices)

    ones_counts = Int[sum(m) for m in matrices]

    check_sizes = expansion_sizes(TARGET_N)
    @printf("[ber] Target N=%d  check_sizes=%s  x_multiplier=%.2f\n",
            TARGET_N, join(check_sizes, ","), X_MULTIPLIER)

    nslots   = Threads.maxthreadid()
    rngs     = [Xoshiro(GLOBAL_SEED + t) for t in 1:nslots]
    buffers  = [Matrix{BigInt}(undef, TARGET_N, TARGET_N) for _ in 1:nslots]
    @printf("[ber] threads=%d\n", Threads.nthreads())

    results = Dict{Symbol, Any}[]
    t_start = time()
    np      = length(P_VALUES)

    for (pi, p) in enumerate(P_VALUES)
        t_p = time()

        p_n, p_4, p_cond, success, total_ext, fail_by_k =
            run_single_p(matrices, ones_counts, p, check_sizes, TARGET_N,
                         X_MULTIPLIER, rngs, buffers)

        push!(results, Dict(
            :p           => p,
            :probability => p_n,
            :p4_sns      => p_4,
            :conditional => p_cond,
            :extensions  => total_ext,
            :successes   => success,
        ))

        elapsed_p = time() - t_p

        @printf("[ber] [%3d/%-3d] p=%.2f  P=%.14f  (cond=%.6f p4=%.6f ext=%d)  %.1fs\n",
                pi, np, Float64(p), Float64(p_n), Float64(p_cond), Float64(p_4), total_ext, elapsed_p)
    end

    elapsed_total = time() - t_start
    @printf("[ber] All %d p‑values done in %.1f s (%.1f min)\n",
            np, elapsed_total, elapsed_total / 60)

    write_json(results)
end

main()
