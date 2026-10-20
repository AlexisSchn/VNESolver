




function set_up_virtual_decomposition(graph::Graph, partition::Vector{Vector{Int}})

    # Setting the decomposition
    v_subgraphs = Subgraph[]
    cut_edges = Edge[]
    nb_appearance = zeros(Int, nv(graph)) # kind of dangerous that its shared, no?? TODO maybe change this...

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


    cut_edge_set = Set(cut_edges)

    for v_nodes in partition

        sort!(v_nodes)
        idx_of_nodes = zeros(Int, nv(graph))
        for (i, v) in enumerate(v_nodes)
            idx_of_nodes[v] = i
            nb_appearance[v] += 1
        end

        subgraph_edges = Edge[]
        cut_edges_src = [Edge[] for _ in vertices(graph)]
        cut_edges_dst = [Edge[] for _ in vertices(graph)]

        for e in edges(graph)
            src_in = idx_of_nodes[src(e)] != 0
            dst_in = idx_of_nodes[dst(e)] != 0
            src_in && dst_in && push!(subgraph_edges, e)
            if e ∈ cut_edge_set
                src_in && push!(cut_edges_src[src(e)], e)
                dst_in && push!(cut_edges_dst[dst(e)], e)
            end
        end
        
        
        push!(v_subgraphs, Subgraph(
            v_nodes,
            subgraph_edges,
            cut_edges_src,
            cut_edges_dst,
            idx_of_nodes,
            Int[],                    # overlapping_nodes
            zeros(Int, nv(graph)),    # idx_overlapping
            nb_appearance,
        ))

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


