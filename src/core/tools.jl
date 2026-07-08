# src/core/tools.jl
# define some useful functions





function get_routing_cost(instance::Instance, edge_routing::Vector{Vector{Int}})
    
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, sdir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    cost = 0

    for (i_e, v_edge) in enumerate(edges(v_g))
        nodes_of_path = edge_routing[i_e]
        for i_node in 1:(length(edge_routing[i_e])-1)
            cost += ve_dem[src(v_edge), dst(v_edge)] * se_cost[nodes_of_path[i_node], nodes_of_path[i_node+1]]
        end
    end

    return cost
end


function get_placement_cost(instance::Instance, node_placement::Vector{Int})
    
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, sdir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    cost = 0

    for v_node in vertices(v_g)
        cost += vn_dem[v_node] * sn_cost[node_placement[v_node]]
    end

    return cost
end


function get_mapping_cost(instance::Instance, mapping::Mapping)
    return get_placement_cost(instance, mapping.node_placement) + get_routing_cost(instance, mapping.edge_routing) 
end