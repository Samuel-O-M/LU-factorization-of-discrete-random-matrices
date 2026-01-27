using Printf
using Base.Threads

const N_BASE = 8
const TARGET_N = 30
const INPUT_FILE = "data/counts_8x8.bin"
const OUTPUT_FILE = "plot_data/upper_bound.json"

struct Signature
    r::Vector{Int}
    c::Vector{Int}
    count::UInt64
    ones::Int
end

function main()
    if !isfile(INPUT_FILE)
        println("Error: File $INPUT_FILE not found.")
        return
    end

    setprecision(BigFloat, 512)

    println("Loading signatures from $INPUT_FILE...")
    signatures = Vector{Signature}()
    
    open(INPUT_FILE, "r") do f
        while !eof(f)
            chunk = read(f, 24)
            if length(chunk) < 24; break; end

            r = Int.(chunk[1:8])
            c = Int.(chunk[9:16])
            count = reinterpret(UInt64, chunk[17:24])[1]
            ones_cnt = sum(r)

            push!(signatures, Signature(r, c, count, ones_cnt))
        end
    end
    println("Loaded $(length(signatures)) signatures.")

    mkpath("plot_data")
    
    p_steps = collect(0.00:0.01:1.00)
    results = Vector{String}(undef, length(p_steps))
    total_steps = length(p_steps)
    completed = Atomic{Int}(0)

    println("Computing bounds for p = 0.00 ... 1.00 using BigFloat (512 bits)...")
    
    Threads.@threads for i in 1:length(p_steps)
        p_val = p_steps[i]
        
        p = BigFloat(p_val)
        q = 1.0 - p
        m = min(p, q)
        M = max(p, q)
        p2q2 = p^2 + q^2

        dim_weights = [(p^k) * (q^(N_BASE - k)) for k in 0:N_BASE]
        sig_weights = [(p^k) * (q^(64 - k)) for k in 0:64]

        coeff_A = Vector{BigFloat}(undef, TARGET_N + 1)
        coeff_B = Vector{BigFloat}(undef, TARGET_N + 1)
        coeff_C = Vector{BigFloat}(undef, TARGET_N + 1)

        for k in (N_BASE + 1):TARGET_N
            k_big = BigFloat(k)
            term_decay = (m > 0) ? m^(k_big - 9) : BigFloat(0.0)
            
            coeff_A[k] = 1.0 - 2 * (q^k_big) + (q^(2 * k_big - 1))
            coeff_B[k] = term_decay * ((q^k_big) - p2q2)
            coeff_C[k] = M^(2 * k_big - 17)
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
            
            for k in (N_BASE + 1):TARGET_N
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
        print("\rProgress: $c / $total_steps")
    end
    
    println()

    open(OUTPUT_FILE, "w") do f
        write(f, "[\n")
        write(f, join(results, ",\n"))
        write(f, "\n]")
    end
    
    println("Done. Data saved to $OUTPUT_FILE")
end

main()