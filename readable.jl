using Glob
using Printf

mkpath("readable_data")

for filename in glob("data/families_*x*.bin")
    base = basename(filename)
    n = 0
    
    try
        n = parse(Int, split(split(base, '_')[2], 'x')[1])
    catch
        continue
    end

    entries = String[]
    open(filename, "r") do f
        while !eof(f)
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

    output_name = joinpath("readable_data", replace(base, ".bin" => ".json"))
    open(output_name, "w") do f
        write(f, "[\n")
        write(f, join(entries, ",\n"))
        write(f, "\n]")
    end
end

n8_file = "data/counts_8x8.bin"
if isfile(n8_file)
    entries = String[]
    open(n8_file, "r") do f
        while !eof(f)
            chunk = read(f, 24)
            if length(chunk) < 24; break; end
            
            r = Int.(chunk[1:8])
            c = Int.(chunk[9:16])
            cnt = reinterpret(UInt64, chunk[17:24])[1]
            
            push!(entries, "  { \"row_sums\": $r, \"col_sums\": $c, \"count\": $cnt }")
        end
    end

    output_name = joinpath("readable_data", replace(basename(n8_file), ".bin" => ".json"))
    open(output_name, "w") do f
        write(f, "[\n")
        write(f, join(entries, ",\n"))
        write(f, "\n]")
    end
end