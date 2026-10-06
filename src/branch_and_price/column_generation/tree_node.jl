

function solve_column_generation_node!(
        model_master::Model, 
        columns, 
        instance::Instance, 
        v_decomposition::OverlappingVirtualDecomposition, 
        pricers, 
        tree_node::TreeNode,
        tree_queue::PriorityQueue{TreeNode, Float64},
        nb_node_created::Int
    )

    time_beginning = time()
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # Set to 0 the columns not possible
    columns_to_unfix=Column[]
    for v_subgraph in v_decomposition.subgraphs
        for v_node in v_subgraph.nodes
            if v_node ∈ keys(tree_node.branching.placement)
                for column in columns[v_subgraph]
                    if column.submapping.node_placement[v_node] != tree_node.branching.placement[v_node]
                        fix(column.variable, 0.; force=true)
                        push!(columns_to_unfix, column)
                    end
                end
            end
            if v_node ∈ keys(tree_node.branching.zoning)
                for column in columns[v_subgraph]
                    if column.submapping.node_placement[v_node] ∉ tree_node.branching.zoning[v_node]
                        fix(column.variable, 0.; force=true)
                        push!(columns_to_unfix, column)
                    end
                end
            end
        end
    end

    optimize!(model_master)

    dual_costs = DualValues(model_master)

    # Column generation

    time_max            = 1000
    nb_iter_max         = 1000
    gap_min             = 0.05
    stabilization_coeff = 0.
    keep_on             = true

    time_overall        = time() - time_beginning
    iter                = 0
    gap                 = 1.
    time_pricer         = 0
    time_rmp            = 0
    rmp_value, lg_bound = 0., tree_node.lower_bound
    nb_total_columns    = 0
    nb_pricer           = 0
    for v_subgraph in v_decomposition.subgraphs
        if length(v_subgraph.nodes) > 1
            nb_pricer += 1
        end
    end


    #println("Starting column generation...")

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

        for v_subgraph in v_decomposition.subgraphs
            if length(v_subgraph.nodes) == 1
                continue
            end

            t=time()
        
            submapping, reduced_cost = solve_greedy_pricer_branchig(instance, v_subgraph, current_dual_costs, shortest_paths.dists, modified_se_cost, tree_node.branching, nb_greedy = 50, time_max = 10)
                
        
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
        #@printf("Iter %-3d; RMP value: %8.3f;   LG bound: %6.3f,   Pricer: %-8s;   New columns %-2d;   Total columns %-4d;   Aver. red. %5.3f;   Time %5.3f;    Gap %2.3f\n", 
        #    iter, rmp_value, lg_bound, "greedy", nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, gap)

        if average_reduced_cost > -0.5 && rmp_value < 10000 && nb_total_columns > 50
            keep_on=false
        end

    end


    #println("MILP TIME!")

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

            submapping, reduced_cost = update_solve_pricer!(pricers[v_subgraph], v_subgraph, dual_costs, instance, tree_node.branching)
            time_pricer += time()-t

            if reduced_cost < -0.0001 && !isnothing(submapping)
                add_column!(model_master, columns, v_subgraph, submapping, instance)
                keep_on = true
                nb_new_cols += 1
                nb_total_columns += 1
            end

            if isinf(reduced_cost)
                # let's run it again to check...?

                new_pricer = set_up_pricer(instance, v_subgraph)
                submapping, reduced_cost = update_solve_pricer!(new_pricer, v_subgraph, dual_costs, instance, tree_node.branching)

                println("Well after all, I got: $reduced_cost")
            end

            total_reduced += reduced_cost
        end

        new_lg_bound = rmp_value + total_reduced
        if new_lg_bound > lg_bound && rmp_value < 10e4
            lg_bound = new_lg_bound
        end

        gap = (rmp_value-lg_bound)/rmp_value
        average_reduced_cost = total_reduced / nb_pricer

        time_overall = time() - time_beginning
        #@printf("Iter %-3d; RMP value: %8.3f;   LG bound: %6.3f,   Pricer: %-8s;   New columns %-2d;   Total columns %-4d;   Aver. red. %5.3f;   Time %5.3f;    Gap %2.3f\n", 
        #    iter, rmp_value, lg_bound, "milp", nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, gap)

    end 

    #println("Finished with $lg_bound after $time_overall with $time_pricer for pricer and $time_rmp for RMP")


    # BRANCHING PART...
    # Has to be done here, i guess :(
    optimize!(model_master)
    if objective_value(model_master) > 10e4
        # Unfixing the stuff
        unique!(columns_to_unfix)
        for column in columns_to_unfix
            unfix(column.variable)
            set_lower_bound(column.variable, 0.0)
            set_upper_bound(column.variable, 1.)
        end
        return lg_bound, Inf
    end

    (v_node, s_node) = first_idea(model_master, instance, v_decomposition, columns)

    # Unfixing the stuff
    unique!(columns_to_unfix)
    for column in columns_to_unfix
        unfix(column.variable)
        set_lower_bound(column.variable, 0.0)
        set_upper_bound(column.variable, 1.)
    end

    if v_node == 0
        println("\n\nNew solution found: $rmp_value but $lg_bound...\n\n\n")
        return lg_bound, rmp_value
    end

    # left
    branching_left = deepcopy(tree_node.branching)
    branching_left.placement[v_node] = s_node
    if v_node ∈ keys(branching_left.zoning)
        delete!(branching_left.zoning, v_node)
    end
    nb_node_created+=1
    node = TreeNode(nb_node_created, tree_node.depth+1, branching_left, lg_bound)
    enqueue!(tree_queue, node, node.lower_bound)

    # right
    branching_right = deepcopy(tree_node.branching)
    if v_node ∈ keys(branching_right.zoning)
        branching_right.zoning[v_node] = setdiff(branching_right.zoning[v_node], [s_node])
    else
        branching_right.zoning[v_node] = setdiff(vertices(s_g), [s_node])
    end
    nb_node_created+=1
    node = TreeNode(nb_node_created, tree_node.depth+1, branching_right, lg_bound)
    enqueue!(tree_queue, node, node.lower_bound)



    return lg_bound, Inf

end


