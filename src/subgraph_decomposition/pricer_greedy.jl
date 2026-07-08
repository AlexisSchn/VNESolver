


# heuristic/tools.jl
# Define some struct and functions that are useful for the heuristics

struct DynamicWeightMatrix{M, C, T} <: AbstractMatrix{T}
    base_dist::M
    capacities::C
    demand::T
end

# Implement the minimum required Interface for Graphs.jl distance matrix
Base.size(d::DynamicWeightMatrix) = size(d.base_dist)
@inline function Base.getindex(d::DynamicWeightMatrix, u::Int, v::Int)
    if d.capacities[u, v] < d.demand
        return Inf
    else
        return d.base_dist[u, v]
    end
end




function shortest_path_routing!(edge_routing, instance::Instance, v_node_placement, v_subgraph::Subgraph, duals::DualValues)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs
    
    # Allocations
    se_cap_copy = copy(se_cap) # TODO : changing?
    nodes_of_path = Vector{Int}()
    sizehint!(nodes_of_path, nv(s_g))
    
    # Tools
    routing_cost    = 0
    virtual_edges   = collect(v_subgraph.edges)

    # Loop
    idx_edges       = Vector(1:length(virtual_edges))
    shuffle!(idx_edges) # routing in a random order. TODO consider the demand/topology

    for i_edge in idx_edges

        if !isempty(edge_routing[i_edge])
            continue
        end

        v_edge = virtual_edges[i_edge]
        demand_curr_edge = ve_dem[src(v_edge), dst(v_edge)]
        idx_src = v_subgraph.idx_of_nodes[src(v_edge)]        
        idx_dst = v_subgraph.idx_of_nodes[dst(v_edge)]
        s_src = v_node_placement[idx_src]
        s_dst = v_node_placement[idx_dst]
         
        current_weights = DynamicWeightMatrix(se_cost, se_cap_copy, demand_curr_edge) # virtual matrix, much faster

        edges_of_path = a_star(s_dir, s_src, s_dst, current_weights)

        if isempty(edges_of_path)
            return Inf  
        end

        empty!(edge_routing[i_edge])
        push!(edge_routing[i_edge], src(edges_of_path[1]))
        cost_of_routing_current_edge = 0
        for edge in edges_of_path
            u, v = src(edge), dst(edge)
            push!(edge_routing[i_edge], v)
            if u < v
                cost_of_routing_current_edge += se_cost[u, v] - duals.edge_capacity[Edge(u, v)]
            else
                cost_of_routing_current_edge += se_cost[u, v] - duals.edge_capacity[Edge(v, u)]
            end
            se_cap_copy[u, v] -= demand_curr_edge
            se_cap_copy[v, u] -= demand_curr_edge # undir version
        end

        routing_cost += cost_of_routing_current_edge * demand_curr_edge
    end

    return routing_cost

end


function complete_partial_placement!(partial_placement::Vector{Int}, instance::Instance, shortest_paths, v_subgraph::Subgraph, duals::DualValues)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # Structures
    nb_nodes = length(v_subgraph.nodes)
    is_placed               = falses(nb_nodes)
    in_frontier             = falses(nb_nodes)
    is_available            = trues(nv(s_g))
    scores                  = zeros(Float64, nv(s_g))
    frontier                = Vector{Int}()
    sizehint!(frontier, nv(v_g))
    placement_v_neighbors   = Vector{Int}()
    sizehint!(placement_v_neighbors, nb_nodes)
    
    # Loop
    placement_cost  = 0

    # Initialization
    for (i_node, v_node) in enumerate(v_subgraph.nodes)
        if !iszero(partial_placement[i_node])
            is_placed[i_node] = true
            curr_placement = partial_placement[i_node]
            is_available[curr_placement] = false
            placement_cost += vn_dem[v_node] * sn_cost[curr_placement] - duals.node_1t1[curr_placement]
            placement_cost -= sum( duals.flow_conservation[curr_placement, v_edge] + duals.flow_departure[curr_placement, v_edge] for v_edge in v_subgraph.cut_edges_with_src[i_node]; init=0.)
            placement_cost += sum( duals.flow_conservation[curr_placement, v_edge] for v_edge in v_subgraph.cut_edges_with_dst[i_node];init=0.)
        end
    end

    if !any(is_placed)
        # TODO tackle that case: place a random node wherever!
    end

    for (i_node, v_node) in enumerate(v_subgraph.nodes)
        if is_placed[i_node]
            for neighbor in neighbors(v_g, v_node)
                if neighbor ∈ v_subgraph.nodes
                    i_neighbor = v_subgraph.idx_of_nodes[neighbor]
                    if !is_placed[i_neighbor] && !in_frontier[i_neighbor]
                        in_frontier[i_neighbor] = true
                        push!(frontier, i_neighbor)
                    end
                end
            end
        end
    end

    while !isempty(frontier)

        shuffle!(frontier)
        i_v_node = pop!(frontier)
        v_node = v_subgraph.nodes[i_v_node]
        curr_demand = vn_dem[v_node]

        empty!(placement_v_neighbors)
        for v_neigh in neighbors(v_g, v_node)
            if v_neigh ∈ v_subgraph.nodes
                i_neighbor = v_subgraph.idx_of_nodes[v_neigh]
                if is_placed[i_neighbor]
                    push!(placement_v_neighbors, partial_placement[i_neighbor])
                end
            end
        end
        
        @. scores = ifelse(is_available & (sn_cap >= curr_demand), 0.0, Inf)
        for p_neigh in placement_v_neighbors
            @views scores .+= shortest_paths.dists[:,p_neigh ]
        end
        selected_node = argmin(scores)

        if isinf(scores[selected_node])
            return Inf
        end

        partial_placement[i_v_node] = selected_node
        placement_cost += sn_cost[selected_node] * curr_demand - duals.node_1t1[selected_node]
        placement_cost -= sum( duals.flow_conservation[selected_node, v_edge] + duals.flow_departure[selected_node, v_edge] for v_edge in v_subgraph.cut_edges_with_src[i_v_node];init=0.)
        placement_cost += sum( duals.flow_conservation[selected_node, v_edge] for v_edge in v_subgraph.cut_edges_with_dst[i_v_node];init=0.)
        is_placed[i_v_node] = true
        is_available[selected_node] = false
        for v_neigh in neighbors(v_g, v_node)
            if v_neigh ∈ v_subgraph.nodes
                i_neighbor = v_subgraph.idx_of_nodes[v_neigh]
                if !is_placed[i_neighbor] && !in_frontier[i_neighbor]
                    in_frontier[i_neighbor] = true
                    push!(frontier, i_neighbor)
                end
            end
        end
    end

    return placement_cost
end
        




function solve_greedy_pricer(instance::Instance, v_subgraph::Subgraph, duals::DualValues; nb_greedy = 100, time_max = 10)

    time_beginning = time()

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # Memory allocation
    placement = Vector{Int}(undef, length(v_subgraph.nodes))
    routing = Vector{Vector{Int}}(undef, length(v_subgraph.edges))
    for i in 1:length(routing)
        routing[i] = Vector{Int}()
        sizehint!(routing[i], nv(s_g))
    end

    best_placement = Vector{Int}(undef, length(v_subgraph.nodes))
    best_routing = Vector{Vector{Int}}(undef, length(v_subgraph.edges))
    for i in 1:length(best_routing)
        best_routing[i] = Vector{Int}()
        sizehint!(best_routing[i], nv(s_g))
    end

    # Ueful tools
    shortest_paths = floyd_warshall_shortest_paths(s_dir, se_cost) # for evaluation of placements
    i_start_node = rand(1:length(v_subgraph.nodes))
    possible_start_s_node = Vector{Int}()
    for s_node in vertices(s_g)
        if sn_cap[s_node] >= vn_dem[v_subgraph.nodes[i_start_node]]
            push!(possible_start_s_node, s_node)
        end
    end
    
    # Loop tools
    best_cost       = Inf
    iter            = 1
    time_overall = time() - time_beginning

    while iter <= nb_greedy && time_overall < time_max

        s_node_start = rand(possible_start_s_node)
        placement .= 0
        for i_edge in 1:length(v_subgraph.edges)
            empty!(routing[i_edge])
        end

        placement[i_start_node] = s_node_start
        placement_cost = complete_partial_placement!(placement, instance, shortest_paths, v_subgraph, duals) 

        if placement_cost < Inf
            routing_cost = shortest_path_routing!(routing, instance, placement, v_subgraph, duals)
        else
            routing_cost = Inf
        end

        total_cost = placement_cost + routing_cost

        if total_cost < best_cost
            best_cost = total_cost
            #println("New best mapping with cost $best_cost")
            best_placement .= placement
            for i_edge in 1:length(v_subgraph.edges)
                empty!(best_routing[i_edge])
                append!(best_routing[i_edge], routing[i_edge])
            end
        end

        iter += 1
        time_overall = time() - time_beginning
    end

    
    if isinf(best_cost)
        return nothing, 0.0
    end

    reduced_cost = best_cost - duals.submapping_selection[v_subgraph]


    #println("Found $reduced_cost with $iter iterations and $(time()-time_beginning) time")
    return Mapping(best_placement, best_routing), reduced_cost
end

