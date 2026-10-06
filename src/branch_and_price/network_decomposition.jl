


struct Subgraph
    nodes::Vector{Int}
    edges::Vector{Edge}
    cut_edges_src::Vector{Vector{Edge}}
    cut_edges_dst::Vector{Vector{Edge}}
    idx_of_nodes::Vector{Int}
    overlapping_nodes::Vector{Int}
    idx_overlapping::Vector{Int}
    nb_appearance_nodes::Vector{Int}
end


struct OverlappingVirtualDecomposition
    subgraphs::Vector{Subgraph}
    cut_edges::Vector{Edge}
    overlapping_nodes::Vector{Int}
    nb_appearance_nodes::Vector{Int}
end

function set_up_virtual_decomposition(graph::Graph, partition::Vector{Vector{Int}})

    # Setting the decomposition
    v_subgraphs = Subgraph[]
    cut_edges = Edge[]
    nb_appearance = zeros(Int, nv(graph))

    for v_edge in edges(graph)
        is_cut_edge = true
        for part in partition
            if src(v_edge) ∈ part && dst(v_edge) ∈ part 
                is_cut_edge = false
            end
        end
        if is_cut_edge
            push!(cut_edges, v_edge)
        end
    end

    for v_nodes in partition
        sort!(v_nodes)
        idx_of_nodes = zeros(Int, nv(graph))
        for (i_node, v_node) in enumerate(v_nodes) 
            idx_of_nodes[v_node] = i_node
            nb_appearance[v_node] += 1
        end


        subgraph_edges = Edge[]
        cut_edges_src = [Edge[] for _ in vertices(graph)]
        cut_edges_dst = [Edge[] for _ in vertices(graph)]

        for v_edge in edges(graph)
            if src(v_edge) ∈ v_nodes && dst(v_edge) ∈ v_nodes
                push!(subgraph_edges, v_edge)
            end
            if src(v_edge) ∈ v_nodes && v_edge ∈ cut_edges
                push!(cut_edges_src[src(v_edge)], v_edge)  
            end
            if dst(v_edge) ∈ v_nodes && v_edge ∈ cut_edges
                push!(cut_edges_dst[dst(v_edge)], v_edge)  
            end
        end
        push!(v_subgraphs, Subgraph(v_nodes, subgraph_edges, cut_edges_src, cut_edges_dst, idx_of_nodes, [], zeros(Int, nv(graph)), nb_appearance))
    end
    
    
    # Overlapping stuff
    overlapping_nodes = Int[]
    for v_node in vertices(graph)
        if nb_appearance[v_node] > 1
            push!(overlapping_nodes, v_node)
        end
    end

    for v_node in overlapping_nodes
        current_i = 1
        for v_subgraph in v_subgraphs
            if v_node in v_subgraph.nodes
                push!(v_subgraph.overlapping_nodes, v_node)
                v_subgraph.idx_overlapping[v_node] = current_i
                current_i += 1
            end
        end
    end
    
    # TODO : check that no edges are in two subgraphs.

    return OverlappingVirtualDecomposition(v_subgraphs, cut_edges, overlapping_nodes, nb_appearance)

end




function set_up_substrate_subgraphs(graph::Graph, partition::Vector{Vector{Int}})

    # Setting the decomposition
    s_subgraphs = Subgraph[]
    cut_edges = Edge[]
    for v_nodes in partition
        sort!(v_nodes)
        idx_of_nodes = zeros(Int, nv(graph))
        for (i_node, v_node) in enumerate(v_nodes) 
            idx_of_nodes[v_node] = i_node
        end

        subgraph_edges = Edge[]
        cut_edges_with_src = [Edge[] for _ in v_nodes]
        cut_edges_with_dst = [Edge[] for _ in v_nodes]

        for v_edge in edges(graph)
            if src(v_edge) ∈ v_nodes
                if dst(v_edge) ∈ v_nodes
                    push!(subgraph_edges, v_edge)
                else
                    push!(cut_edges, v_edge)
                    v_source = findfirst(==(src(v_edge)), v_nodes) 
                    push!(cut_edges_with_src[v_source], v_edge)  
                end
            end
            if dst(v_edge) ∈ v_nodes && src(v_edge) ∉ v_nodes
                v_dst = findfirst(==(dst(v_edge)), v_nodes) 
                push!(cut_edges_with_dst[v_dst], v_edge)  
            end
        end
        push!(s_subgraphs, Subgraph(v_nodes, subgraph_edges, cut_edges_with_src, cut_edges_with_dst, idx_of_nodes, [], [], []))
    end
    
    
    return s_subgraphs

end


function print_stuff_subgraphs(subgraphs::Vector{Subgraph})

    println("There is $(length(subgraphs)) subgraphs:")
    for (i_subgraph, subgraph) in enumerate(subgraphs)
        println("       subgraph $i_subgraph with $(length(subgraph.nodes)) nodes and $(length(subgraph.edges)) edges")
    end
    
end


