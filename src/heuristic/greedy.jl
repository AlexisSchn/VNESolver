# src/heuristic/greedy.jl
# An optimized greedy heuristic, where node placement is achieved depending on distance to previously placed nodes.


function solve_greedy(instance; nb_greedy = 100, time_max = 10)

    time_beginning = time()

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # Memory allocation
    placement = Vector{Int}(undef, nv(v_g))
    routing = Vector{Vector{Int}}(undef, ne(v_g))
    for i in 1:length(routing)
        routing[i] = Vector{Int}()
        sizehint!(routing[i], nv(s_g))
    end

    best_placement = Vector{Int}(undef, nv(v_g))
    best_routing = Vector{Vector{Int}}(undef, ne(v_g))
    for i in 1:length(best_routing)
        best_routing[i] = Vector{Int}()
        sizehint!(best_routing[i], nv(s_g))
    end

    # Ueful tools
    shortest_paths = floyd_warshall_shortest_paths(s_dir, se_cost) # for evaluation of placements
    most_central_v_node = argmin(closeness_centrality(v_g))
    possible_start_s_node = Vector{Int}()
    for s_node in vertices(s_g)
        if sn_cap[s_node] >= vn_dem[most_central_v_node]
            push!(possible_start_s_node, s_node)
        end
    end
    
    # Loop tools
    best_cost       = Inf
    iter            = 1
    time_placement  = 0
    time_routing    = 0
    time_overall = time() - time_beginning


    while iter <= nb_greedy && time_overall < time_max

        s_node_start = rand(possible_start_s_node)
    
        placement .= 0
        for i_edge in 1:ne(v_g)
            empty!(routing[i_edge])
        end

        placement[most_central_v_node] = s_node_start

        time_0 = time()
        placement_cost = complete_partial_placement!(placement, instance, shortest_paths) 
        time_placement += time() - time_0


        time_1 = time()
        if placement_cost < Inf
            routing_cost = shortest_path_routing!(routing, instance, placement)
        else
            routing_cost = Inf
        end
        time_routing += time() - time_1

        total_cost = placement_cost + routing_cost

        if total_cost < best_cost
            best_cost = total_cost
            println("New best mapping with cost $best_cost")
            best_placement .= placement
            for i_edge in 1:ne(v_g)
                empty!(best_routing[i_edge])
                append!(best_routing[i_edge], routing[i_edge])
            end
        end

        iter += 1
        time_overall = time() - time_beginning
    end

    println("Found $best_cost with $iter iterations and $(time()-time_beginning) time, placement $time_placement, routing $time_routing")
    return (mapping=Mapping(best_placement, best_routing), mapping_cost=best_cost)
end



