using Printf
using Base.Threads

const TARGET_N = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 30
const ONES_FILE = "results/ones.json"

struct OnesCount
    ones::Int
    count::UInt128
end

function load_ones_counts(filename::String, n::Int)
    text = read(filename, String)
    section = match(Regex("\"$n\"\\s*:\\s*\\[([^\\]]*)\\]"), text)
    counts = OnesCount[]
    if section !== nothing
        for line in split(section.captures[1], '\n')
            m = match(Regex("\\{\\s*\"k\":\\s*(\\d+),\\s*\"count\":\\s*(\\d+)\\s*\\}"), line)
            if m !== nothing
                push!(counts, OnesCount(parse(Int, m.captures[1]), parse(UInt128, m.captures[2])))
            end
        end
    end
    return counts
end

function compute_bound(ones_counts::Vector{OnesCount}, n_base::Int, target_n::Int, output_file::String)
    setprecision(BigFloat, 512)

    p_steps = collect(0.00:0.01:1.00)
    results = Vector{String}(undef, length(p_steps))
    total_steps = length(p_steps)
    completed = Atomic{Int}(0)

    println("  Computing bounds for p = 0.00 ... 1.00 using BigFloat (512 bits)...")

    Threads.@threads for i in 1:length(p_steps)
        p_val = p_steps[i]

        p = BigFloat(p_val)
        q = 1.0 - p
        m = min(p, q)
        M = max(p, q)
        p2q2 = p^2 + q^2

        sig_weights = [(p^k) * (q^(n_base^2 - k)) for k in 0:(n_base^2)]

        k_values = collect((n_base + 1):target_n)
        in_T = Bool[]
        for k in k_values
            k_big = BigFloat(k)
            lhs = p2q2 * m^(k_big - (n_base + 1))
            rhs = q^k_big * M^(k_big - (n_base + 1)) + (p + (n_base - 1) * m) * M^(2 * k_big - (n_base + 2))
            push!(in_T, lhs >= rhs)
        end

        total_bound = BigFloat(0.0)

        for oc in ones_counts
            w = oc.ones
            w_A = sig_weights[w + 1]
            if w_A == 0.0; continue; end

            group_weight = BigFloat(oc.count) * w_A

            g_w = n_base * exp((BigFloat(w) / n_base) * log(p)) * exp((n_base - BigFloat(w) / n_base) * log(q))

            survival_prob = BigFloat(1.0)

            for (idx, k) in enumerate(k_values)
                if !in_T[idx]; continue; end
                k_big = BigFloat(k)
                term = 2 * q^k_big - q^(2 * k_big - 1) + 2 * (p2q2 * m^(k_big - (n_base + 1)) - q^k_big * M^(k_big - (n_base + 1))) * g_w - M^(2 * k_big - (2 * n_base + 1)) * g_w^2

                if term < 0.0; term = BigFloat(0.0); end
                if term > 1.0; term = BigFloat(1.0); end

                survival_prob *= 1.0 - term

                if survival_prob < 1e-100; break; end
            end

            total_bound += group_weight * survival_prob
        end

        results[i] = @sprintf("  {\"p\": %.2f, \"upper_bound\": %.100e}", Float64(p), total_bound)

        c = atomic_add!(completed, 1) + 1
        print("\r  Progress: $c / $total_steps")
    end

    println()

    open(output_file, "w") do f
        write(f, "[\n")
        write(f, join(results, ",\n"))
        write(f, "\n]")
    end

    println("  Saved to $output_file")
end

function main()
    mkpath("results")

    if !isfile(ONES_FILE)
        println(stderr, "Error: $ONES_FILE not found. Run analyze_ones.jl first.")
        exit(1)
    end

    ones_8 = load_ones_counts(ONES_FILE, 8)
    println("N=8: Loaded $(length(ones_8)) ones counts.")
    compute_bound(ones_8, 8, TARGET_N, "results/weak_upper_bound_8x8.json")
    println("N=8: Done.")

    ones_9 = load_ones_counts(ONES_FILE, 9)
    println("N=9: Loaded $(length(ones_9)) ones counts.")
    compute_bound(ones_9, 9, TARGET_N, "results/weak_upper_bound_9x9.json")
    println("N=9: Done.")
end

main()
