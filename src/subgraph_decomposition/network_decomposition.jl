


struct Subgraph
    nodes::Vector{Int}
    edges::Vector{Edge}
    cut_edges_with_src::Vector{Vector{Edge}}
    cut_edges_with_dst::Vector{Vector{Edge}}
    idx_of_nodes::Vector{Int}
end


struct VirtualDecomposition
    subgraphs::Vector{Subgraph}
    cut_edges::Vector{Edge}
end

function set_up_virtual_decomposition(graph::Graph, partition::Vector{Vector{Int}})

    # Setting the decomposition
    v_subgraphs = Subgraph[]
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
        push!(v_subgraphs, Subgraph(v_nodes, subgraph_edges, cut_edges_with_src, cut_edges_with_dst, idx_of_nodes))
    end
    
    
    return VirtualDecomposition(v_subgraphs, cut_edges)

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
        push!(s_subgraphs, Subgraph(v_nodes, subgraph_edges, cut_edges_with_src, cut_edges_with_dst, idx_of_nodes))
    end
    
    
    return s_subgraphs

end

function print_stuff_subgraphs(subgraphs::Vector{Subgraph})

    println("There is $(length(subgraphs)) subgraphs:")
    for (i_subgraph, subgraph) in enumerate(subgraphs)
        println("       subgraph $i_subgraph with $(length(subgraph.nodes)) nodes and $(length(subgraph.edges)) edges")
    end
    
end



function star_partition(graph::Graph)

    copy_v_network = copy(graph)
    keep_on = true

    real_indices = collect(1:nv(graph)) 
    # Graphs.jl is completly stupid when it comes to removing a node. It swaps it with the last node and then removes it. be careful... 

    nodes_in_no_subgraphs = collect(1:nv(graph))
    v_node_partitionning = Vector{Vector{Int}}()
    possible_centers = collect(1:nv(graph))

    while keep_on

        # Get node with max degree
        nodes_degrees = [degree(copy_v_network, v_node) for v_node in vertices(copy_v_network)]
        node_sorted_degree = sortperm(nodes_degrees, rev=true)
        center = 0
        for node in node_sorted_degree
            if real_indices[node] ∈ possible_centers
                center = node
                break
            end
        end
        if degree(copy_v_network, center) <= 1
            break
        end

        new_part = [ real_indices[center] ]
        for neigh in neighbors(copy_v_network, center)
            #println("Neighbor: $neigh, aka $(real_indices[neigh])")
            if length(new_part) < 8
                push!(new_part, real_indices[neigh])
            end
        end

        for v_node in new_part
            filter!(x -> x != v_node, possible_centers)
        end
    
        sort!(new_part)
        push!(v_node_partitionning, new_part)

        for v_node in new_part
            if v_node ∈ nodes_in_no_subgraphs
                filter!(x -> x != v_node, nodes_in_no_subgraphs)
            end
        end

        while length(neighbors(copy_v_network, center)) > 0
            neigh = neighbors(copy_v_network, center)[1]
            rem_vertex!(copy_v_network, neigh)

            real_indices[neigh] = real_indices[length(real_indices)]
            if center == length(real_indices)
                center = neigh
            end
            deleteat!(real_indices, length(real_indices))


        end

        rem_vertex!(copy_v_network, center)
        real_indices[center] = real_indices[length(real_indices)]
        deleteat!(real_indices, length(real_indices))

    end


    for v_node_left in nodes_in_no_subgraphs
        push!(v_node_partitionning, [v_node_left])
    end

    return v_node_partitionning

end

