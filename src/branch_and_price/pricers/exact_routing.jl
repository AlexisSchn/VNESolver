# ============================================================================
#   Exact routing for a fixed placement
#
#   Routes the virtual edges of a subgraph whose hosts are already chosen, with
#   joint edge capacities, without a MILP. Used by the exact pricer for the leaves
#   where independent shortest paths overload a substrate edge.
#
#   Idea: start with every virtual edge on its shortest path. If a substrate edge
#   is overloaded, at least one of the virtual edges using it must avoid it in any
#   feasible routing. So create one option per user: "this virtual edge may not use
#   that substrate edge", and re-route only that virtual edge. Always continue from
#   the cheapest option (best-first): the first option without overload is optimal.
#
#   Main function: route_exactly
# ============================================================================
 
 
# One option of the search: which substrate edges each virtual edge may not use,
# and the resulting shortest paths.
struct RoutingOption
    forbidden::Vector{Vector{Tuple{Int,Int}}}   # forbidden[e]: substrate edges (u, v) with u < v that virtual edge e may not use
    paths::Vector{Vector{Int}}                  # paths[e]: current path of virtual edge e, as a list of substrate nodes
    path_costs::Vector{Float64}                 # path_costs[e]: cost of paths[e] for one unit of demand
    total::Float64                              # total routing cost: Σ demand × path cost
end
 
 
# Arc weights seen by one virtual edge: Inf on arcs that cannot carry its demand
# and on the substrate edges it is forbidden to use. Both directions of a substrate
# edge share the same capacity, so a forbidden edge is forbidden in both directions.
struct RoutingWeights <: AbstractMatrix{Float64}
    cost::Matrix{Float64}
    cap::Matrix{Int}
    demand::Int
    forbidden::Vector{Tuple{Int,Int}}
end
 
Base.size(w::RoutingWeights) = size(w.cost)
 
function Base.getindex(w::RoutingWeights, u::Int, v::Int)
    if w.cap[u, v] < w.demand
        return Inf
    end
    if minmax(u, v) in w.forbidden
        return Inf
    end
    return w.cost[u, v]
end
 
 
 
"""
    route_exactly(s_dir, cost, cap, dists, commodities, limit; max_options = 1000)
 
Optimal routing of `commodities`, a list of (src host, dst host, demand), one per virtual
edge, with joint capacities `cap` and costs per unit of demand `cost`. `dists[d]` are the
shortest-path distances on edges with capacity ≥ d (used to guide A*).
 
Only routings cheaper than `limit` are of interest. Returns (status, total_cost, paths):
- `:optimal`     : `paths` is an optimal routing, with `total_cost < limit`;
- `:cannot_beat` : no routing costs less than `limit` (this includes "no routing exists");
- `:gave_up`     : more than `max_options` options were needed; use the MILP instead.
"""
function route_exactly(
        s_dir::SimpleDiGraph{Int},
        cost::Matrix{Float64},
        cap::Matrix{Int},
        dists::Vector{Matrix{Float64}},
        commodities::Vector{Tuple{Int,Int,Int}},
        limit::Float64;
        max_options::Int = 1000,
    )
 
    nb_edges = length(commodities)
    demands  = [d for (_, _, d) in commodities]
 
 
    # ============ FIRST OPTION: every virtual edge on its shortest path, nothing forbidden
    no_forbidden = [Tuple{Int,Int}[] for _ in 1:nb_edges]
    paths        = Vector{Vector{Int}}(undef, nb_edges)
    path_costs   = zeros(nb_edges)
 
    for e in 1:nb_edges
        from, to, d = commodities[e]
        route = shortest_route(s_dir, cost, cap, d, no_forbidden[e], from, to, dists[d])
        if isnothing(route)
            return (:cannot_beat, Inf, nothing)         # this virtual edge cannot be routed at all
        end
        paths[e], path_costs[e] = route
    end
 
    total = 0.0
    for e in 1:nb_edges
        total += demands[e] * path_costs[e]
    end
    if total >= limit
        return (:cannot_beat, Inf, nothing)
    end
 
 
    # ============ BEST-FIRST SEARCH
    # `options` stores every option created; the queue holds (total cost, index in `options`),
    # so popping it always gives the cheapest option not explored yet.
    options = [RoutingOption(no_forbidden, paths, path_costs, total)]
    queue   = BinaryMinHeap{Tuple{Float64,Int}}()
    push!(queue, (total, 1))
    seen    = Set{Vector{Tuple{Int,Int,Int}}}()        # forbidden sets already created (avoid duplicates)
 
    while !isempty(queue)
        total, index = pop!(queue)
        if total >= limit
            return (:cannot_beat, Inf, nothing)         # the cheapest option is already too expensive
        end
        option = options[index]
 
        overloaded_edge, users = find_overloaded_edge(option, demands, cap)
        if isnothing(overloaded_edge)
            return (:optimal, total, option.paths)      # no overload and cheapest of all: optimal
        end
 
        if length(options) >= max_options
            return (:gave_up, Inf, nothing)
        end
 
        # One child per user of the overloaded edge: "this virtual edge may not use it".
        # In any feasible routing, at least one of them avoids it, so the children
        # together cover every feasible routing.
        for e in users
            forbidden = [copy(list) for list in option.forbidden]
            push!(forbidden[e], overloaded_edge)
 
            key = forbidden_key(forbidden)
            if key in seen
                continue                                # same option reached by another branch
            end
            push!(seen, key)
 
            # only virtual edge e changes: re-route it, keep the others
            from, to, d = commodities[e]
            route = shortest_route(s_dir, cost, cap, d, forbidden[e], from, to, dists[d])
            if isnothing(route)
                continue                                # e has no path left: dead end
            end
            new_path, new_cost = route
 
            child_total = total - d * option.path_costs[e] + d * new_cost
            if child_total >= limit
                continue                                # cannot beat the limit
            end
 
            child_paths      = copy(option.paths)
            child_paths[e]   = new_path
            child_costs      = copy(option.path_costs)
            child_costs[e]   = new_cost
 
            push!(options, RoutingOption(forbidden, child_paths, child_costs, child_total))
            push!(queue, (child_total, length(options)))
        end
    end
 
    return (:cannot_beat, Inf, nothing)                 # every option was a dead end
end
 
 
 
# ============ HELPERS
 
# Shortest path from `from` to `to` for one virtual edge with demand `demand`, avoiding the
# substrate edges in `forbidden`. Returns (path as a list of nodes, cost per unit), or nothing.
# A* is guided by `dist_d`, the distances without any forbidden edge: forbidding edges can
# only make paths longer, so these distances never overestimate.
function shortest_route(s_dir, cost, cap, demand, forbidden, from, to, dist_d)
    weights    = RoutingWeights(cost, cap, demand, forbidden)
    path_edges = a_star(s_dir, from, to, weights, u -> dist_d[u, to])
    if isempty(path_edges)
        return nothing
    end
 
    path      = [from]
    path_cost = 0.0
    for edge in path_edges
        push!(path, dst(edge))
        path_cost += cost[src(edge), dst(edge)]
    end
    return path, path_cost
end
 
 
# Finds a substrate edge whose load exceeds its capacity. If several are overloaded,
# returns the one with the fewest users (fewer children to create).
# Returns (edge, users) or (nothing, Int[]) if the routing is feasible.
function find_overloaded_edge(option::RoutingOption, demands::Vector{Int}, cap::Matrix{Int})
    load        = Dict{Tuple{Int,Int}, Int}()
    users       = Dict{Tuple{Int,Int}, Vector{Int}}()
    edges_found = Tuple{Int,Int}[]                      # in order of appearance, for reproducibility
 
    for e in eachindex(option.paths)
        path = option.paths[e]
        for i in 1:length(path)-1
            edge = minmax(path[i], path[i+1])
            if !haskey(load, edge)
                load[edge]  = 0
                users[edge] = Int[]
                push!(edges_found, edge)
            end
            load[edge] += demands[e]
            push!(users[edge], e)
        end
    end
 
    best_edge, best_users = nothing, Int[]
    for edge in edges_found
        if load[edge] > cap[edge[1], edge[2]]
            if isnothing(best_edge) || length(users[edge]) < length(best_users)
                best_edge, best_users = edge, users[edge]
            end
        end
    end
    return best_edge, best_users
end
 
 
# A sorted list of (virtual edge, u, v) describing all forbidden substrate edges.
# Two options with the same forbidden edges get the same key.
function forbidden_key(forbidden::Vector{Vector{Tuple{Int,Int}}})
    key = Tuple{Int,Int,Int}[]
    for e in eachindex(forbidden)
        for (u, v) in forbidden[e]
            push!(key, (e, u, v))
        end
    end
    return sort!(key)
end
 
