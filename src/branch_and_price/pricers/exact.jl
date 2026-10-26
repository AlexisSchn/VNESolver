# ============================================================================
#   EXACT COMBINATORIAL PRICER
#
#   Branch-and-bound over node placements, meant for small subgraphs (≤ 5-6 nodes).
#   Same objective as the MILP pricer (update_solve_pricer!):
#         Σ node terms         (node_rc: costs + all node-related duals, Inf if forbidden)
#       + Σ demand × path cost (modified_se_cost: edge costs - capacity duals)
#       - duals.submapping_selection
#
#   Search: virtual nodes are placed one at a time (most constrained first), cheapest
#   hosts first. For a partial placement, the lower bound is:
#         placed part      : exact node terms + shortest paths between placed endpoints
#       + unplaced nodes   : min-cost assignment to DISTINCT free hosts, each host costed
#                            with its node term + shortest paths to its placed neighbours
#       + unplaced edges   : demand × cheapest arc, for edges between two unplaced nodes
#   Shortest paths use dists[d] (only edges with capacity ≥ d) but ignore JOINT capacity,
#   so the bound is valid. At a leaf, edges are routed one by one on residual capacities:
#     - if the routing costs exactly the shortest paths, the leaf is exact;
#     - otherwise (detour or failure) the column is still kept if feasible and good, and the
#       placement is remembered. After the search, the remembered placements that could
#       still beat the best column are routed exactly with a small routing MILP, so the
#       returned min reduced cost is exact.
#
#   Assumptions:
#     - modified_se_cost ≥ 0 (true: capacity duals are ≤ 0);
#     - the MILP pricer also forbids two virtual nodes of the subgraph on the same host;
#     - every virtual edge of the subgraph has a demand ≥ 1.
#   Uses compute_can_host, compute_node_rc, edge_demand and DynamicWeightMatrix
#   from the greedy pricer file.
# ============================================================================
 
 


# ====== TODO
# Implement non-MILP routing; it seems to be the real bottleneck
# The detour stuff?


# Better bound; especially the unplaced edge part sucks right now...
#   [13, 14, 15]: milp -0.4851 (0.132s) | exact -0.4851 (0.093s, 26 nodes) i also had stuff like this;...

Base.@kwdef mutable struct ExactSearch
    # ---- data
    k::Int                                              # number of virtual nodes
    n_sn::Int                                           # number of substrate nodes
    order::Vector{Int}                                  # order[j]: virtual node placed at position j
    node_rc::Matrix{Float64}                            # node_rc[s, v], Inf if forbidden
    forward::Vector{Vector{Tuple{Int,Int,Bool}}}        # forward[i]: (j, demand, j is the edge's src), for later neighbours j
    tail_lb::Vector{Float64}                            # tail_lb[i]: bound for edges whose two ends come after position i
    edges::Vector{Tuple{Int,Int,Int}}                   # (position of src, position of dst, demand), in v_subgraph.edges order
    route_order::Vector{Int}                            # at leaves, edges are routed by decreasing demand
    dists::Vector{Matrix{Float64}}                      # dists[d][s, t]: shortest path using only edges with capacity ≥ d
    modified_se_cost::Matrix{Float64}
    se_cap::Matrix{Int}
    s_dir::SimpleDiGraph{Int}
 
    # ---- search state
    lb::Matrix{Float64}                                 # lb[s, j]: cost of position j on host s, given the placed positions
    saved::Matrix{Float64}                              # lb columns saved before an update, restored on backtrack
    saved_offset::Vector{Int}                           # first column of `saved` used by each level
    perm::Matrix{Int}                                   # candidate hosts of each level, by increasing lb
    used::BitVector                                     # hosts already taken
    place::Vector{Int}                                  # place[j]: host of position j (0 if not placed)
    top_host::Matrix{Int}                               # buffers for the assignment bound
    top_cost::Matrix{Float64}
    picked::Vector{Int}
    suffix_min::Vector{Float64}
    residual::Matrix{Int}                               # residual capacities, for routing at leaves
 
    # ---- results
    best::Float64                                       # best total cost known (starts at the cutoff)
    unresolved::Vector{Tuple{Float64, Vector{Int}}} = Tuple{Float64, Vector{Int}}[]   # (lower bound, placement) of leaves whose routing is not proven optimal
    columns::Vector{Tuple{Float64, Vector{Int}, Vector{Vector{Int}}}} = Tuple{Float64, Vector{Int}, Vector{Vector{Int}}}[]
    nb_nodes::Int = 0
    deadline::Float64
    aborted::Bool = false
end
 
 
 
"""
    solve_exact_pricer(instance, v_subgraph, duals, branching, dists, modified_se_cost; kwargs...)
 
Exact pricing for one (small) subgraph. `dists` and `modified_se_cost` must be computed with
the same `duals`. Returns a named tuple:
- `columns`          : up to `max_columns` (Mapping, reduced cost) pairs, best first, all negative;
- `min_reduced_cost` : valid lower bound on the min reduced cost (0.0 if no negative column exists),
                       to use in the Lagrangian bound;
- `exact`            : false if the time limit was hit (then use the MILP for this subgraph);
- `nb_nodes`         : number of search nodes.
 
`incumbent_cost` (optional): total cost (reduced cost + convexity dual) of a known column,
e.g. the greedy's best, to prune from the start.
"""
function solve_exact_pricer(
        instance::Instance,
        v_subgraph::Subgraph,
        duals::DualValues,
        branching::Branching,
        dists::Vector{Matrix{Float64}},
        modified_se_cost::Matrix{Float64};
        max_columns::Int        = 10,
        time_max::Float64       = 10.,
        incumbent_cost::Float64 = Inf,
        max_routing_milps::Int  = 1000000,
    )
    
    println("\nexact solver for $(v_subgraph.nodes)")
    time_beg = time()
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities
    n_sn = nv(s_g)
 
 
    # ============ NODE TERMS (capacities, branching, all node-related duals)
    can_host = compute_can_host(instance, v_subgraph, branching)
    node_rc  = compute_node_rc(instance, v_subgraph, duals, can_host)
 
 
    # ============ PLACEMENT ORDER AND EDGE STRUCTURE
    order  = placement_order(v_subgraph, v_g, ve_dem)
    k      = length(order)
    pos_of = zeros(Int, nv(v_g))
    pos_of[order] .= 1:k
 
    edges   = [(pos_of[src(e)], pos_of[dst(e)], ve_dem[src(e), dst(e)]) for e in v_subgraph.edges]
    forward = [Tuple{Int,Int,Bool}[] for _ in 1:k]
    for (a, b, d) in edges
        earlier, later = minmax(a, b)
        push!(forward[earlier], (later, d, later == a))
    end
 
    # edges between two unplaced nodes cost at least demand × cheapest usable arc
    cheapest = [cheapest_arc(dists[d]) for d in eachindex(dists)]
    tail_lb  = zeros(k)
    for i in 1:k, (a, b, d) in edges
        if min(a, b) > i
            tail_lb[i] += d * cheapest[d]
        end
    end
 
    saved_offset = zeros(Int, k)
    for i in 2:k
        saved_offset[i] = saved_offset[i-1] + length(forward[i-1])
    end
    nb_saved = saved_offset[k] + length(forward[k])
 
    lb = Matrix{Float64}(undef, n_sn, k)
    for j in 1:k
        lb[:, j] .= @view node_rc[:, order[j]]
    end
 
 
    # ============ SEARCH
    conv        = duals.submapping_selection[v_subgraph]
    cutoff_cost = conv - 1e-6                                   # a column is useful iff its cost is below this
 
    st = ExactSearch(
        k                = k,
        n_sn             = n_sn,
        order            = order,
        node_rc          = node_rc,
        forward          = forward,
        tail_lb          = tail_lb,
        edges            = edges,
        route_order      = sortperm([d for (_, _, d) in edges], rev = true),
        dists            = dists,
        modified_se_cost = modified_se_cost,
        se_cap           = se_cap,
        s_dir            = s_dir,
        lb               = lb,
        saved            = Matrix{Float64}(undef, n_sn, max(nb_saved, 1)),
        saved_offset     = saved_offset,
        perm             = Matrix{Int}(undef, n_sn, k),
        used             = falses(n_sn),
        place            = zeros(Int, k),
        top_host         = zeros(Int, k, k),
        top_cost         = fill(Inf, k, k),
        picked           = zeros(Int, k),
        suffix_min       = zeros(k + 1),
        residual         = similar(se_cap),
        best             = min(cutoff_cost, incumbent_cost),
        deadline         = time() + time_max,
    )
 
    _search!(st, 1, 0.0)
 
    println("Time after exploration: $(time()-time_beg) for $(st.nb_nodes) nodes")

    # ============ LEAVES WHOSE ROUTING WAS NOT PROVEN OPTIMAL: route them exactly
    # Every other placement was either evaluated exactly or pruned with a bound ≥ best,
    # so once these are resolved, st.best is the exact minimum.
    nb_routing_milps = 0
    if !st.aborted && !isempty(st.unresolved)
        sort!(st.unresolved, by = first)
        routing = get_routing_model(v_subgraph, st)         # built once per subgraph, then reused
        set_routing_costs!(routing, st)                     # only the objective changes with the duals
        for (lb_cost, place) in st.unresolved
            lb_cost >= st.best && break                     # sorted: the others cannot improve either
            if nb_routing_milps >= max_routing_milps        # too many conflicts: let the MILP pricer do it
                println("We have too many routings!")
                st.aborted = true
                break
            end
            nb_routing_milps += 1
            routed = _route_exactly(routing, st, place)
            isnothing(routed) && continue                   # joint capacities make this placement infeasible
            routing_cost, paths = routed
            node_cost = sum(st.node_rc[place[j], st.order[j]] for j in 1:st.k)
            total = node_cost + routing_cost
            if total < st.best
                st.best = total
                push!(st.columns, (total, place, paths))
            end
        end
    end

    println("There is $nb_routing_milps MILP todo")
 
    # ============ RESULTS
    min_reduced_cost = st.best < cutoff_cost ? st.best - conv : 0.0
 
    results = Tuple{Mapping, Float64}[]
    for (cost, place, paths) in Iterators.reverse(st.columns)   # improving leaves: the last one is the best
        length(results) >= max_columns && break
        node_placement = zeros(Int, nv(v_g))
        for j in 1:k
            node_placement[order[j]] = place[j]
        end
        push!(results, (Mapping(node_placement, paths), cost - conv))
    end
    
    println("Time overall: $(time() - time_beg)")
    
    return (columns = results, min_reduced_cost = min_reduced_cost, exact = !st.aborted,
            nb_nodes = st.nb_nodes, nb_routing_milps = nb_routing_milps)
end
 
 
 
# ============ PREPARATION HELPERS
 
# Most connected node first, then always the node with the most placed neighbours
# (ties: more neighbours in the subgraph, then larger total edge demand).
function placement_order(v_subgraph::Subgraph, v_g, ve_dem)
    in_sub = falses(nv(v_g))
    in_sub[v_subgraph.nodes] .= true
    is_ordered = falses(nv(v_g))
    order = Int[]
 
    while length(order) < length(v_subgraph.nodes)
        best_v, best_key = 0, (-1, -1, -1)
        for v in v_subgraph.nodes
            is_ordered[v] && continue
            nb_placed  = count(u -> is_ordered[u], neighbors(v_g, v))
            sub_degree = count(u -> in_sub[u], neighbors(v_g, v))
            sub_demand = sum((edge_demand(ve_dem, v, u) for u in neighbors(v_g, v) if in_sub[u]); init = 0)
            key = (nb_placed, sub_degree, sub_demand)
            if key > best_key
                best_v, best_key = v, key
            end
        end
        push!(order, best_v)
        is_ordered[best_v] = true
    end
    return order
end
 
 
# Smallest positive shortest-path distance = cheapest single usable arc (Inf if none).
function cheapest_arc(dist::Matrix{Float64})
    c = Inf
    for x in dist
        if 0 < x < c
            c = x
        end
    end
    return c
end
 
 
 
# ============ SEARCH
 
function _search!(st::ExactSearch, i::Int, cost::Float64)
    st.nb_nodes += 1
    if (st.nb_nodes & 1023) == 0 && time() > st.deadline
        st.aborted = true
    end
    st.aborted && return
 
    col  = @view st.lb[:, i]
    perm = @view st.perm[:, i]
    sortperm!(perm, col)
 
    # Bound on everything after position i, valid whatever host i gets:
    # placing i only adds costs (distances ≥ 0) and removes free hosts.
    rest_floor = st.tail_lb[i]
    for j in i+1:st.k
        rest_floor += _cheapest_free(st, j)
    end
 
    for s in perm
        c = col[s]
        isinf(c) && break                                   # sorted: the remaining hosts are forbidden
        cost + c + rest_floor >= st.best && break           # sorted: no remaining host can do better
        st.used[s] && continue
 
        _place!(st, i, s)
        if i == st.k
            _evaluate_leaf!(st, cost + c)
        else
            limit = st.best - cost - c
            if _rest_lb(st, i, limit) < limit
                _search!(st, i + 1, cost + c)
            end
        end
        _unplace!(st, i, s)
        st.aborted && return
    end
end
 
 
function _place!(st::ExactSearch, i::Int, s::Int)
    st.place[i] = s
    st.used[s]  = true
    off = st.saved_offset[i]
    for (t, (j, d, j_is_src)) in enumerate(st.forward[i])
        col_j = @view st.lb[:, j]
        st.saved[:, off + t] .= col_j
        # the path goes from the host of the edge's src to the host of its dst
        dist = j_is_src ? (@view st.dists[d][:, s]) : (@view st.dists[d][s, :])
        col_j .+= d .* dist
    end
end
 
 
function _unplace!(st::ExactSearch, i::Int, s::Int)
    off = st.saved_offset[i]
    for (t, (j, _, _)) in enumerate(st.forward[i])
        st.lb[:, j] .= @view st.saved[:, off + t]
    end
    st.used[s]  = false
    st.place[i] = 0
end
 
 
 
# ============ BOUNDS
 
function _cheapest_free(st::ExactSearch, j::Int)
    m = Inf
    @inbounds for s in 1:st.n_sn
        if !st.used[s] && st.lb[s, j] < m
            m = st.lb[s, j]
        end
    end
    return m
end
 
 
# Lower bound for positions i+1..k and the edges among them, positions 1..i being placed.
# Returns a value ≥ limit as soon as nothing below `limit` can exist.
function _rest_lb(st::ExactSearch, i::Int, limit::Float64)
    m     = st.k - i
    bound = st.tail_lb[i]
 
    # cheapest free hosts of each remaining position
    conflict = false
    for r in 1:m
        _cheapest_free_hosts!(st, i + r, r, m)
        isinf(st.top_cost[1, r]) && return Inf
        bound += st.top_cost[1, r]
        for r2 in 1:r-1
            conflict |= st.top_host[1, r2] == st.top_host[1, r]
        end
    end
    (!conflict || bound >= limit) && return bound
 
    # several positions want the same host: solve the small assignment exactly
    st.suffix_min[m + 1] = 0.0
    for r in m:-1:1
        st.suffix_min[r] = st.suffix_min[r + 1] + st.top_cost[1, r]
    end
    return st.tail_lb[i] + _assignment!(st, m, 1, 0.0, limit - st.tail_lb[i])
end
 
 
# The m cheapest free hosts of position j, sorted, into column r of the buffers.
function _cheapest_free_hosts!(st::ExactSearch, j::Int, r::Int, m::Int)
    th, tc = st.top_host, st.top_cost
    for q in 1:m
        th[q, r] = 0
        tc[q, r] = Inf
    end
    @inbounds for s in 1:st.n_sn
        st.used[s] && continue
        c = st.lb[s, j]
        c < tc[m, r] || continue
        q = m
        while q > 1 && tc[q - 1, r] > c
            th[q, r] = th[q - 1, r]
            tc[q, r] = tc[q - 1, r]
            q -= 1
        end
        th[q, r] = s
        tc[q, r] = c
    end
end
 
 
# Min-cost assignment of the m remaining positions to distinct hosts. Looking at each
# position's m cheapest free hosts is enough: if a position used a host ranked below m,
# one of its m cheapest would be free, and swapping would not increase the cost.
# Returns `limit` if no assignment costs less than `limit`.
function _assignment!(st::ExactSearch, m::Int, r::Int, acc::Float64, limit::Float64)
    r > m && return acc
    best = limit
    for q in 1:m
        h = st.top_host[q, r]
        h == 0 && break
        c = st.top_cost[q, r]
        acc + c + st.suffix_min[r + 1] >= best && break     # sorted costs
        already = false
        for r2 in 1:r-1
            already |= st.picked[r2] == h
        end
        already && continue
        st.picked[r] = h
        best = min(best, _assignment!(st, m, r + 1, acc + c, best))
    end
    return best
end
 
 
 
# ============ LEAVES
 
# All nodes placed: `lb_cost` = node terms + independent shortest paths. Route the edges on
# residual capacities (largest demands first) to get a feasible column and its true cost.
function _evaluate_leaf!(st::ExactSearch, lb_cost::Float64)
    node_cost = 0.0
    for j in 1:st.k
        node_cost += st.node_rc[st.place[j], st.order[j]]
    end
 
    st.residual .= st.se_cap
    paths        = Vector{Vector{Int}}(undef, length(st.edges))
    routing_cost = 0.0
    for e in st.route_order
        a, b, d      = st.edges[e]
        s_src, s_dst = st.place[a], st.place[b]
        weights      = DynamicWeightMatrix(st.modified_se_cost, st.residual, d)
        dist_d       = st.dists[d]
        path_edges   = a_star(st.s_dir, s_src, s_dst, weights, u -> dist_d[u, s_dst])
 
        if isempty(path_edges)
            # this routing order failed, which does not prove the placement infeasible
            push!(st.unresolved, (lb_cost, copy(st.place)))
            return
        end
 
        path = [s_src]
        for edge in path_edges
            u, v = src(edge), dst(edge)
            push!(path, v)
            routing_cost      += d * st.modified_se_cost[u, v]
            st.residual[u, v] -= d
            st.residual[v, u] -= d
        end
        paths[e] = path
    end
 
    total = node_cost + routing_cost
    if total > lb_cost + 1e-9
        # capacities forced a detour: a cheaper routing may exist for this placement
        push!(st.unresolved, (lb_cost, copy(st.place)))
    end
    if total < st.best
        st.best = total
        push!(st.columns, (total, copy(st.place), paths))
    end
end
 
 
# ============ EXACT ROUTING OF CONFLICTED LEAVES
#
# A small unsplittable multicommodity flow MILP with joint capacities, routing all the
# subgraph's edges for a fixed placement. Its structure only depends on the subgraph and
# the substrate, so it is built ONCE per subgraph and reused:
#   - between calls, only the objective changes (the duals);
#   - between leaves, only the right-hand sides of the endpoints' flow conservation change.
 
mutable struct RoutingModel
    model::Model
    y::Dict{Tuple{Int, Edge{Int}}, VariableRef}                 # (edge index, arc) -> variable
    conservation::Dict{Tuple{Int, Int}, ConstraintRef}          # (edge index, substrate node) -> constraint
    endpoints::Vector{Tuple{Int, Int}}                          # current (src host, dst host) of each edge, (0, 0) if unset
end
 
# One routing model per subgraph, kept for the whole run (keyed by object identity).
const ROUTING_MODELS = IdDict{Any, RoutingModel}()
 
get_routing_model(v_subgraph, st::ExactSearch) = get!(() -> build_routing_model(st), ROUTING_MODELS, v_subgraph)
 
 
function build_routing_model(st::ExactSearch)
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    #set_attribute(model, "CPXPARAM_Threads", 1)
 
    # y[(e, arc)] only for arcs whose capacity can carry the edge's demand
    y = Dict{Tuple{Int, Edge{Int}}, VariableRef}()
    for (e, (_, _, d)) in enumerate(st.edges), arc in edges(st.s_dir)
        if st.se_cap[src(arc), dst(arc)] >= d
            y[(e, arc)] = @variable(model, binary = true)
        end
    end
 
    # flow conservation, right-hand sides set per leaf in _route_exactly
    conservation = Dict{Tuple{Int, Int}, ConstraintRef}()
    for e in eachindex(st.edges), s in 1:st.n_sn
        out_s, in_s = AffExpr(0.0), AffExpr(0.0)
        for t in outneighbors(st.s_dir, s)
            var = get(y, (e, Edge(s, t)), nothing)
            isnothing(var) || add_to_expression!(out_s, var)
        end
        for t in inneighbors(st.s_dir, s)
            var = get(y, (e, Edge(t, s)), nothing)
            isnothing(var) || add_to_expression!(in_s, var)
        end
        conservation[(e, s)] = @constraint(model, out_s - in_s == 0)
    end
 
    # joint capacity of each undirected substrate edge
    for arc in edges(st.s_dir)
        u, v = src(arc), dst(arc)
        u < v || continue
        load = AffExpr(0.0)
        for (e, (_, _, d)) in enumerate(st.edges), a in (Edge(u, v), Edge(v, u))
            var = get(y, (e, a), nothing)
            isnothing(var) || add_to_expression!(load, d, var)
        end
        isempty(load.terms) || @constraint(model, load <= st.se_cap[u, v])
    end
 
    return RoutingModel(model, y, conservation, fill((0, 0), length(st.edges)))
end
 
 
# The objective depends on the duals: set once per pricing call.
function set_routing_costs!(routing::RoutingModel, st::ExactSearch)
    @objective(routing.model, Min,
        sum(st.edges[e][3] * st.modified_se_cost[src(a), dst(a)] * var for ((e, a), var) in routing.y; init = 0.0))
end
 
 
# Routes all edges for `place` (hosts by position). Returns (routing cost, paths) or nothing if infeasible.
function _route_exactly(routing::RoutingModel, st::ExactSearch, place::Vector{Int})
 
    # move each edge's endpoints: reset the old right-hand sides first, then set the new ones
    for (e, (pa, pb, _)) in enumerate(st.edges)
        new_src, new_dst = place[pa], place[pb]
        old_src, old_dst = routing.endpoints[e]
        (old_src, old_dst) == (new_src, new_dst) && continue
        old_src > 0 && set_normalized_rhs(routing.conservation[(e, old_src)], 0.0)
        old_dst > 0 && set_normalized_rhs(routing.conservation[(e, old_dst)], 0.0)
        set_normalized_rhs(routing.conservation[(e, new_src)], 1.0)
        set_normalized_rhs(routing.conservation[(e, new_dst)], -1.0)
        routing.endpoints[e] = (new_src, new_dst)
    end
 
    optimize!(routing.model)
    termination_status(routing.model) == MOI.OPTIMAL || return nothing
 
    # read the paths (optimal flows have no cycles, since costs are positive)
    paths = Vector{Vector{Int}}(undef, length(st.edges))
    for (e, (pa, pb, _)) in enumerate(st.edges)
        s_src, s_dst = place[pa], place[pb]
        path, current = [s_src], s_src
        while current != s_dst
            length(path) > st.n_sn && error("_route_exactly: path does not reach its destination")
            next = 0
            for t in outneighbors(st.s_dir, current)
                var = get(routing.y, (e, Edge(current, t)), nothing)
                if !isnothing(var) && value(var) > 0.5
                    next = t
                    break
                end
            end
            next == 0 && error("_route_exactly: broken path")
            push!(path, next)
            current = next
        end
        paths[e] = path
    end
    return objective_value(routing.model), paths
end
 
