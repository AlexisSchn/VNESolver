# heuristic/tools.jl
# Define some struct and functions that are useful for the heuristics

struct DynamicWeightMatrix{M, C, T} <: AbstractMatrix{T}
    base_dist::M
    capacities::C
    demand::T
end

function DynamicWeightMatrix(base_dist, capacities, demand)
    return DynamicWeightMatrix{typeof(base_dist), typeof(capacities), typeof(demand)}(base_dist, capacities, demand)
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




function shortest_path_routing!(edge_routing, instance, v_node_placement)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs
    
    # Allocations
    se_cap_copy = copy(se_cap) # TODO : changing?
    nodes_of_path = Vector{Int}()
    sizehint!(nodes_of_path, nv(s_g))
    
    # Tools
    routing_cost    = 0
    virtual_edges   = collect(edges(v_g))

    # Init
    for (i_edge, v_edge) in enumerate(virtual_edges)
        if !isempty(edge_routing[i_edge])
            nodes_of_path = edge_routing[i_edge]
            for i_node in 1:length(nodes_of_path)-1
                s_src = nodes_of_path[i_node]
                s_dst = nodes_of_path[i_node+1]
                routing_cost += ve_dem[src(v_edge), dst(v_edge)] * se_cost[s_src, s_dst]
                se_cap_copy[s_src, s_dst] -= ve_dem[src(v_edge), dst(v_edge)]
                se_cap_copy[s_dst, s_src] -= ve_dem[src(v_edge), dst(v_edge)]
            end
        end
    end

    # Loop
    idx_edges       = Vector(1:length(virtual_edges))
    shuffle!(idx_edges) # routing in a random order. TODO consider the demand/topology

    for i_edge in idx_edges

        if !isempty(edge_routing[i_edge])
            continue
        end

        v_edge = virtual_edges[i_edge]
        demand_curr_edge = ve_dem[src(v_edge), dst(v_edge)]
        s_src = v_node_placement[src(v_edge)]
        s_dst = v_node_placement[dst(v_edge)]
         
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
            cost_of_routing_current_edge += se_cost[u, v] 
            se_cap_copy[u, v] -= demand_curr_edge
            se_cap_copy[v, u] -= demand_curr_edge # undir version
        end

        routing_cost += cost_of_routing_current_edge * demand_curr_edge
    end

    return routing_cost

end


function complete_partial_placement!(partial_placement, instance, shortest_paths)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # Structures
    is_placed               = falses(nv(v_g))
    in_frontier             = falses(nv(v_g))
    is_available            = trues(nv(s_g))
    scores                  = zeros(Float64, nv(s_g))
    frontier                = Vector{Int}()
    sizehint!(frontier, nv(v_g))
    placement_v_neighbors   = Vector{Int}()
    sizehint!(placement_v_neighbors, nv(v_g))
    
    # Loop
    placement_cost  = 0

    # Initialization
    for v_node in vertices(v_g)
        if !iszero(partial_placement[v_node])
            is_placed[v_node] = true
            is_available[partial_placement[v_node]] = false
            placement_cost += vn_dem[v_node] * sn_cost[partial_placement[v_node]]
        end
    end

    if !any(is_placed)
        # TODO tackle that case: place a random node wherever!
    end

    for v_node in vertices(v_g)
        if is_placed[v_node]
            for neighbor in neighbors(v_g, v_node)
                if !is_placed[neighbor] && !in_frontier[neighbor]
                    in_frontier[neighbor] = true
                    push!(frontier, neighbor)
                end
            end
        end
    end

    while !isempty(frontier)

        shuffle!(frontier)
        v_node = pop!(frontier)
        curr_demand = vn_dem[v_node]

        empty!(placement_v_neighbors)
        for v_neigh in neighbors(v_g, v_node)
            if is_placed[v_neigh]
                push!(placement_v_neighbors, partial_placement[v_neigh])
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

        partial_placement[v_node] = selected_node
        placement_cost += sn_cost[selected_node] * curr_demand 
        is_placed[v_node] = true
        is_available[selected_node] = false
        for v_neigh in neighbors(v_g, v_node)
            if !is_placed[v_neigh] && !in_frontier[v_neigh]
                in_frontier[v_neigh] = true
                push!(frontier, v_neigh)
            end
        end
    end

    return placement_cost
end
        