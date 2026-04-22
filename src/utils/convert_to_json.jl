# Convert binary data files to readable JSON format
# Usage: julia src/utils/convert_to_json.jl [local|super]
#   local - converts data/local/binary/ files (default)
#   super - converts data/super/binary/ files

using Printf

const MODE = length(ARGS) > 0 ? ARGS[1] : "local"

if MODE == "super"
    const DATA_DIR = "data/super/binary"
    const OUTPUT_DIR = "data/super/readable"
    const FAMILIES_RANGE = 1:7
    const COUNTS_N = 9
else
    const DATA_DIR = "data/local/binary"
    const OUTPUT_DIR = "data/local/readable"
    const FAMILIES_RANGE = 1:7
    const COUNTS_N = 8
end

println("Running in $MODE mode")
println("Input directory: $DATA_DIR")
println("Output directory: $OUTPUT_DIR")

mkpath(OUTPUT_DIR)

for n in FAMILIES_RANGE
    filename = joinpath(DATA_DIR, "families_$(n)x$(n).bin")
    
    if !isfile(filename)
        println("Skipping families_$(n)x$(n).bin (not found)")
        continue
    end
    
    entries = String[]
    open(filename, "r") do f
        # In super mode, families files have an 8-byte MAGIC_FOOTER at the end
        stop_position = MODE == "super" ? stat(filename).size - 8 : stat(filename).size
        
        while position(f) < stop_position
            chunk = read(f, 16)
            if length(chunk) < 16; break; end
            
            rep = reinterpret(UInt64, chunk[1:8])[1]
            count = reinterpret(UInt64, chunk[9:16])[1]

            matrix_str = "[\n"
            for i in 0:n-1
                row_vals = String[]
                for j in 0:n-1
                    bit = (rep >> (i * n + j)) & 1
                    push!(row_vals, string(bit))
                end
                
                line = join(row_vals, ", ")
                if i < n - 1
                    matrix_str *= "    " * line * ",\n"
                else
                    matrix_str *= "    " * line * "\n"
                end
            end
            matrix_str *= "  ]"

            push!(entries, "  { \"matrix\": " * matrix_str * ", \"count\": " * string(count) * " }")
        end
    end

    output_name = joinpath(OUTPUT_DIR, "families_$(n)x$(n).json")
    open(output_name, "w") do f
        write(f, "[\n")
        write(f, join(entries, ",\n"))
        write(f, "\n]")
    end
    println("Converted families_$(n)x$(n).bin -> families_$(n)x$(n).json ($(length(entries)) entries)")
end

counts_file = joinpath(DATA_DIR, "counts_$(COUNTS_N)x$(COUNTS_N).bin")
if isfile(counts_file)
    entries = String[]
    open(counts_file, "r") do f
        while !eof(f)
            if COUNTS_N == 8
                chunk = read(f, 24)
                if length(chunk) < 24; break; end
                r = Int.(chunk[1:8])
                c = Int.(chunk[9:16])
                cnt = reinterpret(UInt64, chunk[17:24])[1]
            else
                # 9x9 format: 9 bytes row_sums + 9 bytes col_sums + 16 bytes count (UInt128)
                chunk = read(f, 34)
                if length(chunk) < 34; break; end
                r = Int.(chunk[1:9])
                c = Int.(chunk[10:18])
                cnt = reinterpret(UInt128, chunk[19:34])[1]
            end
            push!(entries, "  { \"row_sums\": $r, \"col_sums\": $c, \"count\": $cnt }")
        end
    end

    output_name = joinpath(OUTPUT_DIR, "counts_$(COUNTS_N)x$(COUNTS_N).json")
    open(output_name, "w") do f
        write(f, "[\n")
        write(f, join(entries, ",\n"))
        write(f, "\n]")
    end
    println("Converted counts_$(COUNTS_N)x$(COUNTS_N).bin -> counts_$(COUNTS_N)x$(COUNTS_N).json ($(length(entries)) entries)")
else
    println("Skipping counts_$(COUNTS_N)x$(COUNTS_N).bin (not found)")
end

println("Done. All files converted to $OUTPUT_DIR")
