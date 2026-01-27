mkpath("plot_data")

counts = Dict{Int, UInt64}()

open("data/counts_8x8.bin", "r") do f
    while !eof(f)
        chunk = read(f, 24)
        if length(chunk) < 24; break; end
        
        k = sum(chunk[1:8])
        c = reinterpret(UInt64, chunk[17:24])[1]
        
        counts[k] = get(counts, k, 0) + c
    end
end

open("plot_data/ones.json", "w") do f
    write(f, "[\n")
    
    sorted_keys = sort(collect(keys(counts)))
    entries = String[]
    
    for k in sorted_keys
        push!(entries, "  {\"k\": $k, \"count\": $(counts[k])}")
    end
    
    write(f, join(entries, ",\n"))
    write(f, "\n]")
end