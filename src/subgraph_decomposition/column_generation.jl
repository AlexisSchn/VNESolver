


function solve_subgraph_decomposition(instance::Instance; parameters::SubgraphDecompositionParameters=SubgraphDecompositionParameters())

    time_beginning = time()
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # Compute the partition
    #v_partition = [[i] for i in 1:nv(v_g)]
    v_partition = star_partition(v_g)

    println("Partition: $v_partition")
    v_decomposition = set_up_virtual_decomposition(instance.v_network.graph, v_partition)

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_decomposition.subgraphs)
    println("   and $(length(v_decomposition.cut_edges)) cutting edges")



    # RMP
    model_master    = Model(CPLEX.Optimizer)
    set_attribute(model_master, "CPXPARAM_LPMethod", 2)
    set_silent(model_master)
    columns         = set_up_master_problem!(model_master, instance, v_decomposition)

    add_dumb_columns!(model_master, v_decomposition)
    add_single_node_columns!(model_master, columns, instance, v_decomposition)
    add_greedy_columns!(model_master, columns, instance, v_decomposition)

    # Column generation
    pricers = Dict{Subgraph, Model}()
    for v_subgraph in v_decomposition.subgraphs
        pricers[v_subgraph] = set_up_pricer(instance, v_subgraph)
    end
    time_max            = parameters.time_max
    nb_iter_max         = parameters.nb_iter_max
    gap_min             = parameters.gap_min
    stabilization_coeff = parameters.stab_coeff
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

    # ------- Greedy part
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

        for v_subgraph in v_decomposition.subgraphs
            if length(v_subgraph.nodes) == 1
                continue
            end

            t=time()
        
            submapping, reduced_cost = solve_greedy_pricer(instance, v_subgraph, current_dual_costs, shortest_paths.dists, modified_se_cost, nb_greedy = 50, time_max = 10)
                
        
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

        if average_reduced_cost > -0.1 && rmp_value < 10000
            keep_on=false
        end

    end


    println("Switching to exact pricers...")
    keep_on=true
    optimize!(model_master)
    dual_costs = DualValues(model_master)

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
            

        rmp_value = objective_value(model_master)


        for v_subgraph in v_decomposition.subgraphs
            if length(v_subgraph.nodes) == 1
                continue
            end

            t=time()                
                
            submapping, reduced_cost = update_solve_pricer!(pricers[v_subgraph], v_subgraph, dual_costs, instance)
            time_pricer += time()-t

            if reduced_cost < -0.0001 && !isnothing(submapping)
                add_column!(model_master, columns, v_subgraph, submapping, instance)
                keep_on = true
                nb_new_cols += 1
                nb_total_columns += 1
            end

            total_reduced += reduced_cost
        end

        new_lg_bound = rmp_value + total_reduced
        if new_lg_bound > lg_bound
            lg_bound = new_lg_bound
        end

        gap = (rmp_value-lg_bound)/rmp_value
        average_reduced_cost = total_reduced / nb_pricer

        time_overall = time() - time_beginning
        @printf("Iter %-3d; RMP value: %8.3f;   LG bound: %6.3f,   Pricer: %-8s;   New columns %-2d;   Total columns %-4d;   Aver. red. %5.3f;   Time %5.3f;    Gap %2.3f\n", 
            iter, rmp_value, lg_bound, "milp", nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, gap)

    end

    println("Found solution $lg_bound in $time_overall with $time_rmp in RMP and $time_pricer in pricers")

    result = SubgraphDecompositionResult(
        instance.v_network.name,
        instance.s_network.name,
        rmp_value, 
        lg_bound,
        gap,
        nb_total_columns,
        iter,
        time()-time_beginning
    )
    println("Result: $result")
    return result
end





function solve_subgraph_decomposition_better(instance::Instance; parameters::SubgraphDecompositionParameters=SubgraphDecompositionParameters())

    time_beginning = time()
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # Compute the partition
    #v_partition = [[i] for i in 1:nv(v_g)]
    v_partition = star_partition(v_g)

    println("Partition: $v_partition")
    v_decomposition = set_up_virtual_decomposition(instance.v_network.graph, v_partition)

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_decomposition.subgraphs)
    println("   and $(length(v_decomposition.cut_edges)) cutting edges")



    # RMP
    model_master    = Model(CPLEX.Optimizer)
    set_attribute(model_master, "CPXPARAM_LPMethod", 2)
    set_silent(model_master)
    columns         = set_up_master_problem!(model_master, instance, v_decomposition)

    add_dumb_columns!(model_master, v_decomposition)
    add_single_node_columns!(model_master, columns, instance, v_decomposition)
    #add_greedy_columns!(model_master, columns, instance, v_decomposition)

    # Column generation
    pricers = Dict{Subgraph, Model}()
    for v_subgraph in v_decomposition.subgraphs
        pricers[v_subgraph] = set_up_pricer(instance, v_subgraph)
    end
    time_max            = parameters.time_max
    nb_iter_max         = parameters.nb_iter_max
    gap_min             = parameters.gap_min
    stabilization_coeff = parameters.stab_coeff
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
    # I want 
    size_max_v_subgraph = maximum(length(v_subgraph.nodes) for v_subgraph in v_decomposition.subgraphs)
    nb_substrate_subgraphs = floor(Int, nv(s_g) / size_max_v_subgraph)
    part = partition_metis(s_g, nb_substrate_subgraphs)
    clusters = [Int[] for i in 1:nb_substrate_subgraphs]
    for s_node in vertices(s_g)
        push!(clusters[part[s_node]], s_node)
    end

    println("Clusters: $clusters")

    # Now, let's augment it!
    nb_nodes_min = size_max_v_subgraph * 3
    complete_clusters!(clusters, s_g, nb_nodes_min)
    println("Clusters: $clusters")

    print("SUBSTRATE DECOMPO:")
    s_subgraphs = set_up_substrate_subgraphs(s_g, clusters)
    print_stuff_subgraphs(s_subgraphs)

    initialization_greedy!(model_master, columns, instance, v_decomposition, s_subgraphs)


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

        for v_subgraph in v_decomposition.subgraphs
            if length(v_subgraph.nodes) == 1
                continue
            end

            for s_subgraph in s_subgraphs

                t=time()
            
                submapping, reduced_cost = solve_greedy_sub_pricer(instance, v_subgraph, s_subgraph, current_dual_costs, shortest_paths.dists, modified_se_cost, nb_greedy = 50, time_max = 10)
                    
            
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
            iter, rmp_value, lg_bound, "greedy", nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, gap)

        if nb_total_columns > 500 || iter > 25
            keep_on=false
        end

    end


    # ------- Greedy part
    println("Switching to th classical greedy thing!")
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

        for v_subgraph in v_decomposition.subgraphs
            if length(v_subgraph.nodes) == 1
                continue
            end

            t=time()
        
            submapping, reduced_cost = solve_greedy_pricer(instance, v_subgraph, current_dual_costs, shortest_paths.dists, modified_se_cost, nb_greedy = 50, time_max = 10)
                
        
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


    println("Switching to exact pricers...")
    keep_on=true
    optimize!(model_master)
    dual_costs = DualValues(model_master)

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
            

        rmp_value = objective_value(model_master)


        for v_subgraph in v_decomposition.subgraphs
            if length(v_subgraph.nodes) == 1
                continue
            end

            t=time()                
                
            submapping, reduced_cost = update_solve_pricer!(pricers[v_subgraph], v_subgraph, dual_costs, instance)
            time_pricer += time()-t

            if reduced_cost < -0.0001 && !isnothing(submapping)
                add_column!(model_master, columns, v_subgraph, submapping, instance)
                keep_on = true
                nb_new_cols += 1
                nb_total_columns += 1
            end

            total_reduced += reduced_cost
        end

        new_lg_bound = rmp_value + total_reduced
        if new_lg_bound > lg_bound
            lg_bound = new_lg_bound
        end

        gap = (rmp_value-lg_bound)/rmp_value
        average_reduced_cost = total_reduced / nb_pricer

        time_overall = time() - time_beginning
        @printf("Iter %-3d; RMP value: %8.3f;   LG bound: %6.3f,   Pricer: %-8s;   New columns %-2d;   Total columns %-4d;   Aver. red. %5.3f;   Time %5.3f;    Gap %2.3f\n", 
            iter, rmp_value, lg_bound, "milp", nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, gap)

    end

    println("Found solution $lg_bound in $time_overall with $time_rmp in RMP and $time_pricer in pricers")

    result = SubgraphDecompositionResult(
        instance.v_network.name,
        instance.s_network.name,
        rmp_value, 
        lg_bound,
        gap,
        nb_total_columns,
        iter,
        time()-time_beginning
    )
    println("Result: $result")
    return result
end



function add_greedy_columns!(model_master::Model, columns, instance::Instance, v_decomposition::VirtualDecomposition)
    
    
    empty_duals = zero_duals(model_master)
    shortest_paths = floyd_warshall_shortest_paths(instance.s_network.directed_graph, instance.s_network.edge_costs)

    for v_subgraph in v_decomposition.subgraphs
        if length(v_subgraph.nodes) <= 2
            continue
        end
        for i in 1:50
            submapping, reduced_cost = solve_greedy_pricer(instance, v_subgraph, empty_duals, shortest_paths.dists, instance.s_network.edge_costs, nb_greedy = 3, time_max = 10)
            if isnothing(submapping)
                continue
            end
            #println("Found a nice mapping $reduced_cost")
            #=
            if column_already_there(columns, submapping, v_subgraph)
                println("already there!")
                continue
            end
            =#
            add_column!(model_master, columns, v_subgraph, submapping, instance)
        end
    end
end



function complete_clusters!(clusters, s_g, nb_nodes_min)

    for cluster in clusters
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



function initialization_greedy!(model_master::Model, columns, instance::Instance, v_decomposition::VirtualDecomposition, s_subgraphs)
    
    
    empty_duals = zero_duals(model_master)
    shortest_paths = floyd_warshall_shortest_paths(instance.s_network.directed_graph, instance.s_network.edge_costs)

    for v_subgraph in v_decomposition.subgraphs
        if length(v_subgraph.nodes) <= 2
            continue
        end

        for s_subgraph in s_subgraphs

            t=time()
        
            submapping, reduced_cost = solve_greedy_sub_pricer(instance, v_subgraph, s_subgraph, empty_duals, shortest_paths.dists, instance.s_network.edge_costs, nb_greedy = 50, time_max = 10)

            if isnothing(submapping)
                continue
            end
            #println("Found a nice mapping $reduced_cost")
            #=
            if column_already_there(columns, submapping, v_subgraph)
                println("already there!")
                continue
            end
            =#
            add_column!(model_master, columns, v_subgraph, submapping, instance)
        end
    end
end
