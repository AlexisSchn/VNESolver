

 # OK LETS GO

# ==== ROUTING STUFF


function shortest_path_routing!(edge_routing, instance::Instance, v_node_placement::Vector{Int}, v_subgraph::Subgraph, modified_se_cost::Matrix, dists)

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
        s_src = v_node_placement[src(v_edge)]
        s_dst = v_node_placement[dst(v_edge)]
         
        current_weights = DynamicWeightMatrix(modified_se_cost, se_cap_copy, demand_curr_edge) # virtual matrix, much faster

        # edges_of_path = a_star(s_dir, s_src, s_dst, current_weights)
        edges_of_path = a_star(s_dir, s_src, s_dst, current_weights, u -> dists[demand_curr_edge][u, s_dst])

        if isempty(edges_of_path)
            return Inf  
        end

        empty!(edge_routing[i_edge])
        push!(edge_routing[i_edge], src(edges_of_path[1]))
        cost_of_routing_current_edge = 0
        for edge in edges_of_path
            u, v = src(edge), dst(edge)
            push!(edge_routing[i_edge], v)
            cost_of_routing_current_edge += modified_se_cost[u, v] 
            se_cap_copy[u, v] -= demand_curr_edge
            se_cap_copy[v, u] -= demand_curr_edge # undir version
        end

        routing_cost += cost_of_routing_current_edge * demand_curr_edge
    end

    return routing_cost

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


function initialize_placement!(placement, v_subgraph, branching, can_host)
    placement .= 0

    has_mandatory = false
    for v in v_subgraph.nodes
        if haskey(branching.placement_mandatory, v)
            placement[v] = branching.placement_mandatory[v]
            has_mandatory = true
        end
    end
    has_mandatory && return true

    # start from the most constrained virtual node, random tie-break
    nb_hosts(v) = count(@view can_host[:, v])
    fewest = minimum(nb_hosts, v_subgraph.nodes)
    fewest == 0 && return false # if no node can host anything

    start = rand([v for v in v_subgraph.nodes if nb_hosts(v) == fewest])
    placement[start] = rand(findall(@view can_host[:, start]))
    return true
end


edge_demand(ve_dem, u, v) = u < v ? ve_dem[u, v] : ve_dem[v, u]

# Most placed neighbours first, random tie-break. Returns an index into `frontier`.
function pick_most_connected(frontier, nb_placed_nbrs)
    best = maximum(v -> nb_placed_nbrs[v], frontier)
    return rand([i for i in eachindex(frontier) if nb_placed_nbrs[frontier[i]] == best])
end

# Random choice among the `rcl_size` best hosts (rcl_size = 1: best host, random tie-break).
# I HAVE NO CLUE WHAT THJIS IS DOING...
function pick_host(scores, rcl_size)
    best = minimum(scores)
    isinf(best) && return 0
    if rcl_size <= 1
        return rand(findall(x -> x <= best + 1e-9, scores))
    end
    k = min(rcl_size, count(isfinite, scores))
    return rand(partialsortperm(scores, 1:k))
end



function complete_partial_placement!(
        placement::Vector{Int},
        instance::Instance,
        v_subgraph::Subgraph,
        node_rc::Matrix{Float64},
        dists::Vector{Matrix{Float64}};
        rcl_size::Int = 1,
    )
    v_g, ve_dem = instance.v_network.graph, instance.v_network.edge_demands
    n_vn, n_sn = nv(v_g), size(node_rc, 1)

    is_placed      = falses(n_vn)
    in_frontier    = falses(n_vn)
    nb_placed_nbrs = zeros(Int, n_vn)
    is_available   = trues(n_sn)
    scores         = Vector{Float64}(undef, n_sn)
    frontier       = Int[]

    function add_to_frontier!(v)
        for u in neighbors(v_g, v)
            (in_subgraph(v_subgraph, u) && !is_placed[u]) || continue
            nb_placed_nbrs[u] += 1
            if !in_frontier[u]
                in_frontier[u] = true
                push!(frontier, u)
            end
        end
    end

    # 1. Register the nodes already placed, then build the frontier
    #    (two passes, so two adjacent pre-placed nodes don't add each other)
    for v in v_subgraph.nodes
        if placement[v] != 0
            is_placed[v] = true
            is_available[placement[v]] = false
        end
    end
    any(v -> is_placed[v], v_subgraph.nodes) || return Inf
    for v in v_subgraph.nodes
        is_placed[v] && add_to_frontier!(v)
    end

    # 2. Grow the placement
    while !isempty(frontier)
        i = pick_most_connected(frontier, nb_placed_nbrs)
        v = frontier[i]
        frontier[i] = frontier[end]
        pop!(frontier)

        # node term (Inf where forbidden or taken) ...
        @views @. scores = ifelse(is_available, node_rc[:, v], Inf)
        # ... + demand-weighted distance to each placed neighbour
        # (assumes symmetric costs: dists[d][s, p] == dists[d][p, s])
        for u in neighbors(v_g, v)
            if in_subgraph(v_subgraph, u) && is_placed[u]
                d = edge_demand(ve_dem, v, u)
                @views scores .+= d .* dists[d][:, placement[u]]
            end
        end

        s = pick_host(scores, rcl_size)
        s == 0 && return Inf

        placement[v] = s
        is_placed[v] = true
        is_available[s] = false
        add_to_frontier!(v)
    end

    # 3. A disconnected subgraph would leave nodes unplaced
    all(v -> is_placed[v], v_subgraph.nodes) || return Inf

    return sum(node_rc[placement[v], v] for v in v_subgraph.nodes)
end











# GREEDY PRICER
# The claude way...
function solve_greedy_pricer_claude(
        instance::Instance, 
        v_subgraph::Subgraph, 
        s_zone::Vector{Int},
        duals::DualValues, 
        dists::Vector{Matrix{Float64}}, 
        modified_se_cost::Matrix, 
        branching::Branching; 
        nb_greedy = 100, 
        time_max = 10
    )
    
    v_g, vn_dem = instance.v_network.graph, instance.v_network.node_demands
    s_g, sn_cap = instance.s_network.graph, instance.s_network.node_capacities

    time_beginning = time()



    # Memory allocation
    placement = Vector{Int}(undef, nv(v_g))
    routing = Vector{Vector{Int}}(undef, length(v_subgraph.edges))
    for i in 1:length(routing)
        routing[i] = Vector{Int}()
        sizehint!(routing[i], nv(s_g))
    end

    best_placement = Vector{Int}(undef, nv(v_g))
    best_routing = Vector{Vector{Int}}(undef, length(v_subgraph.edges))
    for i in 1:length(best_routing)
        best_routing[i] = Vector{Int}()
        sizehint!(best_routing[i], nv(s_g))
    end

    can_host = compute_can_host(instance, v_subgraph, branching; s_zone = s_zone)
    node_rc  = compute_node_rc(instance, v_subgraph, duals, can_host)

    # Loop tools
    best_cost       = Inf
    iter            = 1
    time_overall = time() - time_beginning

    while iter <= nb_greedy && time_overall < time_max
        
        placement .= 0
        for i_edge in 1:length(v_subgraph.edges)
            empty!(routing[i_edge])
        end


        initialize_placement!(placement, v_subgraph, branching, can_host)

        placement_cost = complete_partial_placement!(placement, instance, v_subgraph, node_rc, dists; rcl_size = 2) 

        if placement_cost < Inf
            routing_cost = shortest_path_routing!(routing, instance, placement, v_subgraph, modified_se_cost, dists)
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
