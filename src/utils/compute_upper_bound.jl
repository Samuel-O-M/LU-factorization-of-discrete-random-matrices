using Printf
using Base.Threads

const TARGET_N = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 30

struct Signature
    r::Vector{Int}
    c::Vector{Int}
    count::UInt128
    ones::Int
end

function load_signatures_8(filename::String)
    sigs = Vector{Signature}()
    open(filename, "r") do f
        while !eof(f)
            chunk = read(f, 24)
            if length(chunk) < 24; break; end
            r = Int.(chunk[1:8])
            c = Int.(chunk[9:16])
            count = UInt128(reinterpret(UInt64, chunk[17:24])[1])
            push!(sigs, Signature(r, c, count, sum(r)))
        end
    end
    return sigs
end

function load_signatures_9(filename::String)
    sigs = Vector{Signature}()
    open(filename, "r") do f
        while !eof(f)
            chunk = read(f, 34)
            if length(chunk) < 34; break; end
            r = Int.(chunk[1:9])
            c = Int.(chunk[10:18])
            count = reinterpret(UInt128, chunk[19:34])[1]
            push!(sigs, Signature(r, c, count, sum(r)))
        end
    end
    return sigs
end

function compute_bound(signatures::Vector{Signature}, n_base::Int, target_n::Int, output_file::String)
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

        dim_weights = [(p^k) * (q^(n_base - k)) for k in 0:n_base]
        sig_weights = [(p^k) * (q^(n_base^2 - k)) for k in 0:(n_base^2)]

        coeff_A = Vector{BigFloat}(undef, target_n + 1)
        coeff_B = Vector{BigFloat}(undef, target_n + 1)
        coeff_C = Vector{BigFloat}(undef, target_n + 1)

        for k in (n_base + 1):target_n
            k_big = BigFloat(k)
            term_decay = (m > 0) ? m^(k_big - (n_base + 1)) : BigFloat(0.0)

            coeff_A[k] = 1.0 - 2 * (q^k_big) + (q^(2 * k_big - 1))
            coeff_B[k] = term_decay * ((q^k_big) - p2q2)
            coeff_C[k] = M^(2 * k_big - (2 * n_base + 1))
        end

        total_bound = BigFloat(0.0)

        for sig in signatures
            w_A = sig_weights[sig.ones + 1]
            if w_A == 0.0; continue; end

            group_weight = BigFloat(sig.count) * w_A

            sum_P_row = sum(dim_weights[val + 1] for val in sig.r)
            sum_P_col = sum(dim_weights[val + 1] for val in sig.c)

            sum_P_total = sum_P_row + sum_P_col
            prod_P_total = sum_P_row * sum_P_col

            survival_prob = BigFloat(1.0)

            for k in (n_base + 1):target_n
                term = coeff_A[k] + (coeff_B[k] * sum_P_total) + (coeff_C[k] * prod_P_total)

                if term < 0.0; term = BigFloat(0.0); end
                if term > 1.0; term = BigFloat(1.0); end

                survival_prob *= term

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

    file_8 = "data/local/binary/counts_8x8.bin"
    file_9 = "data/super/binary/counts_9x9.bin"

    if isfile(file_8)
        println("N=8: Loading signatures from $file_8...")
        sigs_8 = load_signatures_8(file_8)
        println("N=8: Loaded $(length(sigs_8)) signatures.")
        compute_bound(sigs_8, 8, TARGET_N, "results/upper_bound_8x8.json")
        println("N=8: Done.")
    else
        println("N=8: SKIPPED ($file_8 not found)")
    end

    if isfile(file_9)
        println("N=9: Loading signatures from $file_9...")
        sigs_9 = load_signatures_9(file_9)
        println("N=9: Loaded $(length(sigs_9)) signatures.")
        compute_bound(sigs_9, 9, TARGET_N, "results/upper_bound_9x9.json")
        println("N=9: Done.")
    else
        println("N=9: SKIPPED ($file_9 not found)")
    end
end

main()
