


function edge_decomposition(graph::Graph)

    partition = Vector{Vector{Int}}()
    for e in edges(graph)
        push!(partition, [src(e), dst(e)])
    end

    return partition
end