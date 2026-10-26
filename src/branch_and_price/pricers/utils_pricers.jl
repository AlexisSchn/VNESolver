# ============================================================================
#   Util functions for (so far) the homemade VNE algos (exact and greedy
#   Routing: bastract matrix useful for the shortest path routing
#   Placement: can_host and node_rc are useful   
#
#   @notes this is work-in-progress! 
#   @author Alexis Schneider
# ============================================================================




# ==== ROUTING STUFF
struct DynamicWeightMatrix <: AbstractMatrix{Float64}
    base_dist::Matrix{Float64}    
    capacities::Matrix{Int}                   
    demand::Int                                  
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



# ==== PLACEMENT STUFF
# can_host[s, v] == true  ⇔  substrate node s may host virtual node v
function compute_can_host(instance, v_subgraph, branching; s_zone = vertices(instance.s_network.graph))
    vn_dem = instance.v_network.node_demands
    sn_cap = instance.s_network.node_capacities
    n_sn = nv(instance.s_network.graph)

    in_zone = falses(n_sn)
    in_zone[s_zone] .= true

    can_host = falses(n_sn, nv(instance.v_network.graph))
    for v in v_subgraph.nodes
        # If it is mandatory placement, we don't check for the capacity; 
        # it's suppose to be good either way!
        if haskey(branching.placement_mandatory, v) 
            can_host[branching.placement_mandatory[v], v] = true
        elseif haskey(branching.placement_allowed, v)
            for s in branching.placement_allowed[v]
                can_host[s, v] = in_zone[s] && sn_cap[s] >= vn_dem[v]
            end
        else
            can_host[:, v] .= in_zone .& (sn_cap .>= vn_dem[v])
        end
    end
    return can_host
end


# node_rc[s, v] = reduced cost of placing v on s  (Inf if forbidden)
function compute_node_rc(instance, v_subgraph, duals, can_host)
    vn_dem  = instance.v_network.node_demands
    sn_cost = instance.s_network.node_costs
    node_rc = fill(Inf, size(can_host))

    for v in v_subgraph.nodes
        k = v_subgraph.nb_appearance_nodes[v]
        for s in axes(can_host, 1)
            can_host[s, v] || continue
            c = (vn_dem[v] * sn_cost[s] - duals.node_1t1[s]) / k
            for e in v_subgraph.cut_edges_src[v]
                c -= (duals.flow_conservation[s, e] + duals.flow_departure[s, e]) / k
            end
            for e in v_subgraph.cut_edges_dst[v]
                c += duals.flow_conservation[s, e] / k
            end
            if v ∈ v_subgraph.overlapping_nodes
                idx = v_subgraph.idx_overlapping[v]
                nxt = idx < k ? idx + 1 : 1
                c += duals.overlapping[s, v, idx] - duals.overlapping[s, v, nxt]
            end
            node_rc[s, v] = c
        end
    end
    return node_rc
end



edge_demand(ve_dem, u, v) = u < v ? ve_dem[u, v] : ve_dem[v, u]

