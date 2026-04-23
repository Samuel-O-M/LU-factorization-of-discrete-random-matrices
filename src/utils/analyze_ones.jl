# Analyze ones distribution from counts data
# Usage: julia src/utils/analyze_ones.jl [local|super]
#   local - uses data/local/binary/counts_8x8.bin (default)
#   super - uses data/super/binary/counts_9x9.bin

const MODE = length(ARGS) > 0 ? ARGS[1] : "local"

if MODE == "super"
    const INPUT_FILE = "data/super/binary/counts_9x9.bin"
    const OUTPUT_DIR = "results/super"
    const OUTPUT_FILE = joinpath(OUTPUT_DIR, "ones_9x9.json")
    const N = 9
else
    const INPUT_FILE = "data/local/binary/counts_8x8.bin"
    const OUTPUT_DIR = "results/local"
    const OUTPUT_FILE = joinpath(OUTPUT_DIR, "ones_8x8.json")
    const N = 8
end

println("Running in $MODE mode")
println("Input: $INPUT_FILE")
println("Output: $OUTPUT_FILE")

mkpath(OUTPUT_DIR)

counts = Dict{Int, UInt128}()

open(INPUT_FILE, "r") do f
    while !eof(f)
        if N == 8
            chunk = read(f, 24)
            if length(chunk) < 24; break; end
            k = sum(chunk[1:8])
            c = reinterpret(UInt64, chunk[17:24])[1]
        else
            # 9x9 format: 9 bytes row_sums + 9 bytes col_sums + 16 bytes count (UInt128)
            chunk = read(f, 34)
            if length(chunk) < 34; break; end
            k = sum(chunk[1:9])
            c = reinterpret(UInt128, chunk[19:34])[1]
        end
        counts[k] = get(counts, k, UInt128(0)) + c
    end
end

open(OUTPUT_FILE, "w") do f
    write(f, "[\n")
    sorted_keys = sort(collect(keys(counts)))
    entries = String[]
    for k in sorted_keys
        push!(entries, "  {\"k\": $k, \"count\": $(counts[k])}")
    end
    write(f, join(entries, ",\n"))
    write(f, "\n]")
end

println("Done. Results saved to $OUTPUT_FILE")
