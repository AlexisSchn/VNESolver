

function solve_column_generation_root_node!(
        model_master::Model,
        columns, instance::Instance, 
        v_decomposition::OverlappingVirtualDecomposition, 
        pricers
    )

    time_beginning = time()
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    empty_branching = Branching()
    max_dem = maximum(ve_dem)


    optimize!(model_master)

    dual_costs = DualValues(model_master)


    # Column generation

    time_max            = 1000
    nb_iter_max         = 1000
    gap_min             = 0.005
    stabilization_coeff = 0.0
    keep_on             = true

    time_overall        = time() - time_beginning
    iter                = 0
    gap                 = 1.
    time_pricer         = 0
    time_rmp            = 0
    rmp_value, lg_bound = 0.,0.
    nb_total_columns    = 0
    nb_pricer           = 0
    for v_subgraph in v_decomposition.subgraphs
        if length(v_subgraph.nodes) > 1
            nb_pricer += 1
        end
    end


    

    # ------- Substrate SN part
    size_max_v_subgraph = maximum(length(v_subgraph.nodes) for v_subgraph in v_decomposition.subgraphs)
    nb_substrate_subgraphs = floor(Int, nv(s_g) / size_max_v_subgraph)
    part = partition_metis(s_g, nb_substrate_subgraphs)
    clusters = [Int[] for i in 1:nb_substrate_subgraphs]
    for s_node in vertices(s_g)
        push!(clusters[part[s_node]], s_node)
    end

    #println("Clusters: $clusters")

    # Now, let's augment it!
    nb_nodes_min = size_max_v_subgraph * 3
    complete_clusters!(clusters, s_g, nb_nodes_min)



    while keep_on && time_overall < time_max && iter < nb_iter_max && gap > gap_min
        
        iter += 1        
        keep_on         = false
        total_reduced   = 0
        nb_new_cols     = 0
        
        t=time()

        optimize!(model_master)
        time_rmp += time()-t

        current_dual_costs = DualValues(model_master)

        rmp_value = objective_value(model_master)

        modified_se_cost = zeros(Float64, ne(s_g), ne(s_g))
        for edge in edges(s_g)
            u, v = src(edge), dst(edge)
            # Subtract the duals to reflect the true reduced cost of the path
            modified_se_cost[u, v] = se_cost[u, v] - current_dual_costs.edge_capacity[Edge(u, v)]
            modified_se_cost[v, u] = se_cost[v, u] - current_dual_costs.edge_capacity[Edge(u, v)]
        end

        dists = [floyd_warshall_shortest_paths(s_dir, DynamicWeightMatrix(modified_se_cost, se_cap, d)).dists
                for d in 1:max_dem]

        for v_subgraph in v_decomposition.subgraphs
            if length(v_subgraph.nodes) == 1
                continue
            end

            for s_zone in clusters

                t=time()
                
                submapping, reduced_cost = solve_greedy_pricer_claude(
                    instance, 
                    v_subgraph, 
                    s_zone, 
                    current_dual_costs, 
                    dists, 
                    modified_se_cost,
                    empty_branching,
                    nb_greedy = 50, 
                    time_max = 10
                )
                    
                time_pricer += time()-t

                if reduced_cost < -0.0001 && !isnothing(submapping)
                    add_column!(model_master, columns, v_subgraph, submapping, instance)
                    keep_on = true
                    nb_new_cols += 1
                    nb_total_columns += 1
                end
                total_reduced += reduced_cost
            end
        end

        gap = (rmp_value-lg_bound)/rmp_value
        average_reduced_cost = total_reduced / nb_pricer*nb_substrate_subgraphs

        time_overall = time() - time_beginning
        @printf("Iter %-3d; RMP value: %8.3f;   LG bound: %6.3f,   Pricer: %-8s;   New columns %-2d;   Total columns %-4d;   Aver. red. %5.3f;   Time %5.3f;    Gap %2.3f\n", 
            iter, rmp_value, lg_bound, "greedy subpricers", nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, gap)

        if nb_total_columns > 1000 || iter > 100 || nb_new_cols < 3
            keep_on=false
        end

    end

    # ------- Greedy part
    keep_on=true
    while keep_on && time_overall < time_max && iter < nb_iter_max && gap > gap_min
        
        iter += 1        
        keep_on         = false
        total_reduced   = 0
        nb_new_cols     = 0
        
        t=time()

        optimize!(model_master)
        time_rmp += time()-t

        current_dual_costs = DualValues(model_master)

        rmp_value = objective_value(model_master)

        modified_se_cost = zeros(Float64, ne(s_g), ne(s_g))
        for edge in edges(s_g)
            u, v = src(edge), dst(edge)
            # Subtract the duals to reflect the true reduced cost of the path
            modified_se_cost[u, v] = se_cost[u, v] - current_dual_costs.edge_capacity[Edge(u, v)]
            modified_se_cost[v, u] = se_cost[v, u] - current_dual_costs.edge_capacity[Edge(u, v)]
        end
        shortest_paths = floyd_warshall_shortest_paths(s_dir, modified_se_cost)

        dists = [floyd_warshall_shortest_paths(s_dir, DynamicWeightMatrix(modified_se_cost, se_cap, d)).dists
                for d in 1:max_dem]


        for v_subgraph in v_decomposition.subgraphs
            if length(v_subgraph.nodes) == 1
                continue
            end

            t=time()
        
            submapping, reduced_cost = solve_greedy_pricer_claude(
                    instance, 
                    v_subgraph, 
                    collect(1:nv(s_g)), 
                    current_dual_costs, 
                    dists, 
                    modified_se_cost,
                    empty_branching,
                    nb_greedy = 100, 
                    time_max = 10
                )
                
        
            time_pricer += time()-t

            if reduced_cost < -0.0001 && !isnothing(submapping)
                add_column!(model_master, columns, v_subgraph, submapping, instance)
                keep_on = true
                nb_new_cols += 1
                nb_total_columns += 1
            end

            total_reduced += reduced_cost
        end

        gap = (rmp_value-lg_bound)/rmp_value
        average_reduced_cost = total_reduced / nb_pricer

        time_overall = time() - time_beginning
        @printf("Iter %-3d; RMP value: %8.3f;   LG bound: %6.3f,   Pricer: %-8s;   New columns %-2d;   Total columns %-4d;   Aver. red. %5.3f;   Time %5.3f;    Gap %2.3f\n", 
            iter, rmp_value, lg_bound, "greedy", nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, gap)

        if average_reduced_cost > -0.5 && rmp_value < 10000 && nb_total_columns > 1000
            keep_on=false
        end

    end



    # --- MILP

    keep_on=true
    optimize!(model_master)
    dual_costs = DualValues(model_master)
    rmp_value = objective_value(model_master)

    while keep_on && time_overall < time_max && iter < nb_iter_max && gap > gap_min

        iter += 1        
        keep_on         = false
        total_reduced   = 0
        nb_new_cols     = 0
        
        t=time()

        optimize!(model_master)
        time_rmp += time()-t

        current_dual_costs = DualValues(model_master)
        stabilize_duals!(dual_costs, current_dual_costs, stabilization_coeff)
        current_rmp_value = objective_value(model_master)
        rmp_value = stabilization_coeff * rmp_value + (1-stabilization_coeff) * current_rmp_value


        modified_se_cost = zeros(Float64, ne(s_g), ne(s_g))
        for edge in edges(s_g)
            u, v = src(edge), dst(edge)
            # Subtract the duals to reflect the true reduced cost of the path
            modified_se_cost[u, v] = se_cost[u, v] - current_dual_costs.edge_capacity[Edge(u, v)]
            modified_se_cost[v, u] = se_cost[v, u] - current_dual_costs.edge_capacity[Edge(u, v)]
        end
        dists = [floyd_warshall_shortest_paths(s_dir, DynamicWeightMatrix(modified_se_cost, se_cap, d)).dists
                for d in 1:max_dem]


        t=time()                

        for v_subgraph in v_decomposition.subgraphs
            if length(v_subgraph.nodes) == 1
                continue
            end

            if length(v_subgraph.nodes) <= 10
                result = solve_exact_pricer(instance, v_subgraph, dual_costs, empty_branching, dists, modified_se_cost)
                if result.exact
                    for (submapping, rc) in result.columns
                        add_column!(model_master, columns, v_subgraph, submapping, instance)
                        nb_new_cols += 1 
                        nb_total_columns += 1
                    end
                    keep_on |= !isempty(result.columns)
                    total_reduced += result.min_reduced_cost
                    continue
                end
            else
                submapping, reduced_cost = update_solve_pricer!(pricers[v_subgraph], v_subgraph, dual_costs, instance, empty_branching)
                if reduced_cost < -0.0001 && !isnothing(submapping)
                    add_column!(model_master, columns, v_subgraph, submapping, instance)
                    keep_on = true
                    nb_new_cols += 1
                    nb_total_columns += 1
                end
            end


            # t_milp  = @elapsed (_, rc_milp) = update_solve_pricer!(pricers[v_subgraph], v_subgraph, dual_costs, instance, empty_branching)
            # t_exact = @elapsed result = solve_exact_pricer(instance, v_subgraph, dual_costs, empty_branching, dists, modified_se_cost)
            # @printf("  %s: milp %.4f (%.3fs) | exact %.4f (%.3fs, %d nodes)\n",
            #     v_subgraph.nodes[1:min(5, end)], rc_milp, t_milp, result.min_reduced_cost, t_exact, result.nb_nodes)



            total_reduced += reduced_cost
        end
            
        time_pricer += time()-t

        new_lg_bound = rmp_value + total_reduced
        if new_lg_bound > lg_bound && rmp_value < 10e4
            lg_bound = new_lg_bound
        end


        gap = (current_rmp_value-lg_bound)/rmp_value
        average_reduced_cost = total_reduced / nb_pricer

        time_overall = time() - time_beginning
        @printf("Iter %-3d; RMP value: %8.3f;   LG bound: %6.3f,   Pricer: %-8s;   New columns %-2d;   Total columns %-4d;   Aver. red. %5.3f;   Time %5.3f;   Time exact %5.3f;   Gap %2.3f\n", 
            iter, rmp_value, lg_bound, "milp", nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, time_pricer, gap)

        if gap<0.15 && stabilization_coeff<0.01
            println("activating stabilization")
            stabilization_coeff = 0.9
        end

    end 

    println("Finished with $lg_bound after $time_overall with $time_pricer for pricer and $time_rmp for RMP")

    return lg_bound, dual_costs

end








function complete_clusters!(clusters, s_g, nb_nodes_min)

    for cluster in clusters
        if length(cluster) > 0
            while length(cluster) < nb_nodes_min
                s_neighbors = Int[]
                for s_node in cluster
                    for s_neigh in neighbors(s_g, s_node)
                        if s_neigh ∉ cluster
                            push!(s_neighbors, s_neigh)
                        end
                    end
                end
                push!(cluster, rand(s_neighbors))
            end
        end
    end

end




function add_edge_columns!(model_master, columns, instance::Instance, v_decomposition::OverlappingVirtualDecomposition)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    nb_cols_added=0

    for v_subgraph in v_decomposition.subgraphs
        if length(v_subgraph.nodes) == 2
            v_edge = v_subgraph.edges[1]
            for s_edge in edges(s_g)
                if se_cap[src(s_edge), dst(s_edge)] >= ve_dem[src(v_edge), dst(v_edge)]
                    if sn_cap[src(s_edge)] >= vn_dem[src(v_edge)] && sn_cap[dst(s_edge)] >= vn_dem[dst(v_edge)]
                        placement = zeros(Int, nv(v_g))
                        placement[src(v_edge)] = src(s_edge)
                        placement[dst(v_edge)] = dst(s_edge)
                        routing = Vector{Vector{Int}}()
                        push!(routing, [src(s_edge), dst(s_edge)])
                        m = Mapping(placement, routing)
                        add_column!(model_master, columns, v_subgraph, m, instance)
                        nb_cols_added+=1
                    end
                    if sn_cap[dst(s_edge)] >= vn_dem[src(v_edge)] && sn_cap[src(s_edge)] >= vn_dem[dst(v_edge)]
                        placement = zeros(Int, nv(v_g))
                        placement[src(v_edge)] = dst(s_edge)
                        placement[dst(v_edge)] = src(s_edge)
                        routing = Vector{Vector{Int}}()
                        push!(routing, [dst(s_edge), src(s_edge)])
                        m = Mapping(placement, routing)
                        add_column!(model_master, columns, v_subgraph, m, instance)
                        nb_cols_added+=1
                    end
                end
            end
        end
    end

    println("$nb_cols_added columns added during init")
end


