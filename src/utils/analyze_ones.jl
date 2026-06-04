const OUTPUT_FILE = "results/ones.json"
const MAGIC_FOOTER = 0xCAFEBABECAFEBABE

mkpath("results")

all_results = Dict{Int, Dict{Int, UInt128}}()

function read_families_ones(filename::String, n::Int)
    counts = Dict{Int, UInt128}()
    fs = stat(filename).size
    data_size = fs
    open(filename, "r") do io
        seek(io, fs - sizeof(UInt64))
        if read(io, UInt64) == MAGIC_FOOTER
            data_size = fs - sizeof(UInt64)
        end
        seek(io, 0)
        num_records = div(data_size, 16)
        mask = (UInt64(1) << (n * n)) - 1
        for _ in 1:num_records
            rep = read(io, UInt64)
            cnt = read(io, UInt64)
            k = Int(count_ones(rep & mask))
            counts[k] = get(counts, k, UInt128(0)) + UInt128(cnt)
        end
    end
    return counts
end

for n in 1:9
    counts = Dict{Int, UInt128}()
    used_file = nothing

    if n <= 7
        local_file = "data/local/binary/families_$(n)x$(n).bin"
        super_file = "data/super/binary/families_$(n)x$(n).bin"
        if isfile(local_file)
            counts = read_families_ones(local_file, n)
            used_file = local_file
        elseif isfile(super_file)
            counts = read_families_ones(super_file, n)
            used_file = super_file
        end
    elseif n == 8
        filename = "data/local/binary/counts_8x8.bin"
        if isfile(filename)
            open(filename, "r") do io
                while !eof(io)
                    chunk = read(io, 24)
                    if length(chunk) < 24; break; end
                    k = sum(chunk[1:8])
                    c = reinterpret(UInt64, chunk[17:24])[1]
                    counts[k] = get(counts, k, UInt128(0)) + UInt128(c)
                end
            end
            used_file = filename
        end
    else
        filename = "data/super/binary/counts_9x9.bin"
        if isfile(filename)
            open(filename, "r") do io
                while !eof(io)
                    chunk = read(io, 34)
                    if length(chunk) < 34; break; end
                    k = sum(chunk[1:9])
                    c = reinterpret(UInt128, chunk[19:34])[1]
                    counts[k] = get(counts, k, UInt128(0)) + c
                end
            end
            used_file = filename
        end
    end

    if used_file !== nothing
        total = sum(values(counts); init=UInt128(0))
        println("N=$n done ($used_file) Total=$total")
    else
        println("N=$n SKIPPED (no data file found)")
    end
    all_results[n] = counts
end

open(OUTPUT_FILE, "w") do f
    write(f, "{\n")
    n_entries = String[]
    for n in 1:9
        counts = all_results[n]
        sorted_keys = sort(collect(keys(counts)))
        k_entries = String[]
        for k in sorted_keys
            push!(k_entries, "      {\"k\": $k, \"count\": $(counts[k])}")
        end
        push!(n_entries, "  \"$n\": [\n" * join(k_entries, ",\n") * "\n  ]")
    end
    write(f, join(n_entries, ",\n"))
    write(f, "\n}\n")
end

println("Done. Results saved to $OUTPUT_FILE")
