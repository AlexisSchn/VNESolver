# src/heuristic/local_search.jl
# iteratively removes a few adjacent nodes in the graph, replace them, and reroute the adjacent edges

function local_search(instance::Instance, initial_mapping::Mapping; 
        nb_local_search::Int = 500, 
        time_max::Float64 = 10., 
        nb_nodes_removed::Int = 10
    )


    time_beginning = time()

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # Memory allocation
    placement = copy(initial_mapping.node_placement)
    routing = Vector{Vector{Int}}(undef, ne(v_g))
    best_placement = copy(placement)
    best_routing = Vector{Vector{Int}}(undef, ne(v_g))
    for i in 1:length(routing)
        routing[i] = copy(initial_mapping.edge_routing[i])
        best_routing[i] = copy(initial_mapping.edge_routing[i])
    end


    # Ueful tools
    shortest_paths = floyd_warshall_shortest_paths(s_dir, se_cost) # for evaluation of placements
    nb_nodes_removed = minimum([nb_nodes_removed, nv(v_g)])
    v_nodes_deleted = zeros(Int, nb_nodes_removed)
    boundary_nodes  = Vector{Int}()
    sizehint!(boundary_nodes, nv(v_g))

    # Loop
    best_cost       = get_placement_cost(instance, placement) + get_routing_cost(instance, routing)
    step            = 1
    time_overall    = time() - time_beginning
    
    println("Before Local Search, we have solution of $best_cost !")
    
    while step <= nb_local_search && time_overall < time_max

        placement .= best_placement
        
        some_v_node = rand(1:nv(v_g))
        v_nodes_deleted .= 0
        empty!(boundary_nodes)
        v_nodes_deleted[1] = some_v_node
        placement[some_v_node] = 0

        for v_neigh in neighbors(v_g, some_v_node)
            push!(boundary_nodes, v_neigh)
        end

        for i_node in 2:nb_nodes_removed
            shuffle!(boundary_nodes)
            next_node = pop!(boundary_nodes)
            placement[next_node] = 0
            v_nodes_deleted[i_node] = next_node
            for neighbor in neighbors(v_g, next_node)
                if neighbor ∉ v_nodes_deleted && neighbor ∉ boundary_nodes
                    push!(boundary_nodes, neighbor)
                end
            end
        end
        placement_cost = complete_partial_placement!(placement, instance, shortest_paths) 

        if !isinf(placement_cost) 

            for (i_edge, v_edge) in enumerate(edges(v_g))
                empty!(routing[i_edge])
                if src(v_edge) ∉ v_nodes_deleted && dst(v_edge) ∉ v_nodes_deleted
                    for i_node in best_routing[i_edge]
                        push!(routing[i_edge], i_node)
                    end
                end
            end

            routing_cost = shortest_path_routing!(routing, instance, placement)

            current_cost = placement_cost + routing_cost
            if current_cost < best_cost
                println("New best! $current_cost, at iter $step")
                best_cost = current_cost
                best_placement .= placement
                for i_edge in 1:ne(v_g)
                    empty!(best_routing[i_edge])
                    append!(best_routing[i_edge], routing[i_edge])
                end
            end
        end
        step += 1
        time_overall = time() - time_beginning
    end

    println("Find the solution of $best_cost in $(time()-time_beginning) and $step iters")
    return (mapping = Mapping(best_placement, best_routing),
            mapping_cost = best_cost)
end