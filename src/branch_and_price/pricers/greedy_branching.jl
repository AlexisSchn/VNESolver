




function complete_partial_placement!(partial_placement::Vector{Int}, instance::Instance, dists::Matrix, v_subgraph::Subgraph, duals::DualValues, branching::Branching)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # Structures
    nb_nodes = nv(v_g)
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
    for v_node in v_subgraph.nodes
        if !iszero(partial_placement[v_node])
            is_placed[v_node] = true
            curr_placement = partial_placement[v_node]
            is_available[curr_placement] = false
        end
    end

    if !any(is_placed)
        # TODO tackle that case: place a random node wherever!
    end

    for v_node in v_subgraph.nodes
        if is_placed[v_node]
            for neighbor in neighbors(v_g, v_node)
                if neighbor ∈ v_subgraph.nodes
                    if !is_placed[neighbor] && !in_frontier[neighbor]
                        in_frontier[neighbor] = true
                        push!(frontier, neighbor)
                    end
                end
            end
        end
    end
    while !isempty(frontier)

        shuffle!(frontier)
        v_node = pop!(frontier)
        curr_demand = vn_dem[v_node]

        empty!(placement_v_neighbors)
        v_neighbors = Int[]

        for v_neigh in neighbors(v_g, v_node)
            if v_neigh ∈ v_subgraph.nodes
                if is_placed[v_neigh]
                    push!(placement_v_neighbors, partial_placement[v_neigh])
                    push!(v_neighbors, v_neigh)
                end
            end
        end
        
        if v_node ∈ keys(branching.zoning)
            is_allowed = falses(nv(s_g))
            for s_node in branching.zoning[v_node]
                is_allowed[s_node] = true
            end
        else
            is_allowed = trues(nv(s_g))
        end

        @. scores = ifelse(is_available & (sn_cap >= curr_demand) & is_allowed, 0.0, Inf)

        for p_neigh in placement_v_neighbors
            #@views scores .+= shortest_paths.dists[:,p_neigh ]
            #v_neigh = v_neighbors[i_neigh]
            #dem_edge = ve_dem[v_node, v_neigh]
            #if v_node > v_neighbors[i_neigh]
            #    dem_edge = ve_dem[v_neigh, v_node]
            #end
            @views scores .+= dists[:,p_neigh] # * dem_edge .+ (sn_cost-duals.node_1t1) * curr_demand
        end
        selected_node = argmin(scores)

        if isinf(scores[selected_node])
            return Inf
        end

        partial_placement[v_node] = selected_node
        is_placed[v_node] = true
        is_available[selected_node] = false
        for v_neigh in neighbors(v_g, v_node)
            if v_neigh ∈ v_subgraph.nodes
                if !is_placed[v_neigh] && !in_frontier[v_neigh]
                    in_frontier[v_neigh] = true
                    push!(frontier, v_neigh)
                end
            end
        end
    end

    placement_cost = 0
    for v_node in v_subgraph.nodes
        curr_placement = partial_placement[v_node]
        placement_cost += (vn_dem[v_node] * sn_cost[curr_placement] - duals.node_1t1[curr_placement])/v_subgraph.nb_appearance_nodes[v_node]
        placement_cost -= sum( duals.flow_conservation[curr_placement, v_edge]/v_subgraph.nb_appearance_nodes[v_node] + duals.flow_departure[curr_placement, v_edge]/v_subgraph.nb_appearance_nodes[v_node] for v_edge in v_subgraph.cut_edges_src[v_node]; init=0.)
        placement_cost += sum( duals.flow_conservation[curr_placement, v_edge]/v_subgraph.nb_appearance_nodes[v_node] for v_edge in v_subgraph.cut_edges_dst[v_node];init=0.)
        if v_node ∈ v_subgraph.overlapping_nodes
            idx = v_subgraph.idx_overlapping[v_node]
            placement_cost += duals.overlapping[curr_placement, v_node, idx]
            if idx < v_subgraph.nb_appearance_nodes[v_node]
                placement_cost -= duals.overlapping[curr_placement, v_node, idx+1]
            else
                placement_cost -= duals.overlapping[curr_placement, v_node, 1]
            end
        end
    end

    return placement_cost
end
        




function solve_greedy_pricer_branchig(instance::Instance, v_subgraph::Subgraph, duals::DualValues, current_dists::Matrix, modified_se_cost::Matrix, branching::Branching; nb_greedy = 100, time_max = 10)
    time_beginning = time()

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

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

    
    # Loop tools
    best_cost       = Inf
    iter            = 1
    time_overall = time() - time_beginning

    while iter <= nb_greedy && time_overall < time_max
        
        placement .= 0
        still_init = true

        for (v_node, s_node) in branching.placement
            if v_node ∈ v_subgraph.nodes
                placement[v_node] = s_node
                still_init = false
            end
        end
        if still_init
            v_node_zoning = Int[]
            for v_node in keys(branching.zoning)
                if v_node ∈ v_subgraph.nodes
                    push!(v_node_zoning, v_node)
                end
            end
            if !isempty(v_node_zoning)
                start_node = rand(keys(branching.zoning))
                possible_start_s_node = Vector{Int}()
                for s_node in branching.zoning[start_node]
                    if sn_cap[s_node] >= vn_dem[start_node]
                        push!(possible_start_s_node, s_node)
                    end
                end
                s_node_start = rand(possible_start_s_node)
                placement[start_node] = s_node_start
            else
                start_node = rand(v_subgraph.nodes)
                possible_start_s_node = Vector{Int}()
                for s_node in vertices(s_g)
                    if sn_cap[s_node] >= vn_dem[start_node]
                        push!(possible_start_s_node, s_node)
                    end
                end
                s_node_start = rand(possible_start_s_node)
                placement[start_node] = s_node_start

            end
        end
        for i_edge in 1:length(v_subgraph.edges)
            empty!(routing[i_edge])
        end

        placement_cost = complete_partial_placement!(placement, instance, current_dists, v_subgraph, duals, branching) 

        if placement_cost < Inf
            routing_cost = shortest_path_routing!(routing, instance, placement, v_subgraph, modified_se_cost)
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

