

module GraphPartitioning

using Graphs
using Libdl

function graph_to_csr(g::Graph)
    n = nv(g)
    xadj = Vector{Int32}()
    adjncy = Vector{Int32}()

    push!(xadj, 0)
    for i in 1:n
        for neighbor in neighbors(g, i)
            push!(adjncy, Int32(neighbor - 1))
        end
        push!(xadj, Int32(length(adjncy)))
    end
    return xadj, adjncy
end

include("kahip_wrapper.jl")
include("metis_wrapper.jl")


export partition_kahip, partition_metis

end # module