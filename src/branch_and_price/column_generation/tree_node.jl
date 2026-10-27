

function solve_column_generation_node!(
        bnp_state::BranchAndPriceState,
        tree_node::TreeNode
    )

    verbose = 1
    time_beginning = time()
    
    (; instance, model_master, columns, pricers, v_decomposition, best_solution) = bnp_state

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs
    
    max_dem = maximum(ve_dem)


    # ============ PRE TREATMENT : fixing variables depending on the branching decisions

    # columns depending on x' branching
    variable_fixed=VariableRef[]
    for v_subgraph in v_decomposition.subgraphs
        for v_node in v_subgraph.nodes
            if v_node ∈ keys(tree_node.branching.placement_mandatory)
                for column in columns[v_subgraph]
                    if column.submapping.node_placement[v_node] != tree_node.branching.placement_mandatory[v_node]
                        fix(column.variable, 0.; force=true)
                        push!(variable_fixed, column.variable)
                    end
                end
            end
            if v_node ∈ keys(tree_node.branching.placement_allowed)
                for column in columns[v_subgraph]
                    if column.submapping.node_placement[v_node] ∉ tree_node.branching.placement_allowed[v_node]
                        fix(column.variable, 0.; force=true)
                        push!(variable_fixed, column.variable)
                    end
                end
            end
        end
    end
    # y variables for branching
    for (v_edge, s_edges) in tree_node.branching.routing_cut_edges_mandatory
        for s_edge in s_edges
            var = model_master[:y][v_edge, s_edge]
            fix(var, 1.0; force = true)
            push!(variable_fixed, var)
        end
    end
    for v_edge in keys(tree_node.branching.routing_cut_edges_forbidden)
        for s_edge in tree_node.branching.routing_cut_edges_forbidden[v_edge]
            var = model_master[:y][v_edge, s_edge]
            fix(var, 0.; force=true)
            push!(variable_fixed, var)
        end
    end



    # ==== AND ALSO: TRYING SOME FIXIN IN THE ALGO?? get the reduced cost for each variable!

    # =============== COLUMN GENERATION

    optimize!(model_master)

    dual_costs = DualValues(model_master)


    time_max            = 1000
    nb_iter_max         = 1000
    gap_min             = 0.01
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



    # greedy pricers
    keep_on=false
    while keep_on && time_overall < time_max && iter < nb_iter_max
        
        iter += 1        
        keep_on         = false
        total_reduced   = 0
        nb_new_cols     = 0
        
        t=time()

        optimize!(model_master)
        time_rmp += time()-t

        current_dual_costs = DualValues(model_master)



        rmp_value = objective_value(model_master)

        modified_se_cost = zeros(Float64, nv(s_g), nv(s_g))
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

            t=time()         
            submapping, reduced_cost = solve_greedy_pricer_claude(
                instance, 
                v_subgraph, 
                collect(1:nv(s_g)), 
                current_dual_costs, 
                dists, 
                modified_se_cost,
                tree_node.branching,
                nb_greedy = 250, 
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

        average_reduced_cost = total_reduced / nb_pricer

        time_overall = time() - time_beginning
        # @printf("Iter %-3d; RMP value: %8.3f;   LG bound: %6.3f,   Pricer: %-8s;   New columns %-2d;   Total columns %-4d;   Aver. red. %5.3f;   Time %5.3f;    Gap %2.3f\n", 
        #     iter, rmp_value, lg_bound, "greedy", nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, gap)

        if average_reduced_cost > -0.5 && rmp_value < 10000 && nb_total_columns > 50
            keep_on=false
        end

    end


    # exact pricers

    keep_on=true
    optimize!(model_master)
    dual_costs = DualValues(model_master)
    rmp_value = objective_value(model_master)
    heuristic_mode = true
    while keep_on && time_overall < time_max && iter < nb_iter_max && gap > gap_min

        iter += 1        
        keep_on         = false
        total_reduced   = 0
        nb_new_cols     = 0
        
        time_rmp += @elapsed optimize!(model_master)

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

            # t_milp  = @elapsed (_, rc_milp) = update_solve_pricer!(pricers[v_subgraph], v_subgraph, dual_costs, instance, empty_branching)
            # t_exact = @elapsed result = solve_exact_pricer(instance, v_subgraph, dual_costs, empty_branching, dists, modified_se_cost)
            # @printf("  %s: milp %.4f (%.3fs) | exact %.4f (%.3fs, %d nodes)\n",
            #     v_subgraph.nodes[1:min(10, end)], rc_milp, t_milp, result.min_reduced_cost, t_exact, result.nb_nodes)

            use_milp = length(v_subgraph.nodes) > 10

            if !use_milp
                result = solve_exact_pricer(instance, v_subgraph, dual_costs, tree_node.branching,
                                            dists, modified_se_cost; heuristic_mode = heuristic_mode)
                if result.exact
                    for (submapping, rc) in result.columns
                        add_column!(model_master, columns, v_subgraph, submapping, instance)
                        nb_new_cols += 1
                        nb_total_columns += 1
                    end
                    keep_on |= !isempty(result.columns)
                    total_reduced += result.min_reduced_cost
                else
                    @warn "Exact pricer gave up on $(v_subgraph.nodes), falling back to the MILP"
                    use_milp = true
                end
            end

            if use_milp
                submapping, reduced_cost = update_solve_pricer!(pricers[v_subgraph], v_subgraph, dual_costs,
                                                                instance, tree_node.branching)
                if reduced_cost < -0.0001 && !isnothing(submapping)
                    add_column!(model_master, columns, v_subgraph, submapping, instance)
                    keep_on = true
                    nb_new_cols += 1
                    nb_total_columns += 1
                end
                total_reduced += reduced_cost
            end
        end
            
        time_pricer             += time()-t
        new_lg_bound            = rmp_value + total_reduced
        if !heuristic_mode && new_lg_bound > lg_bound && rmp_value < 10e4
            lg_bound = new_lg_bound
        end
        gap                     = (current_rmp_value-lg_bound)/rmp_value
        average_reduced_cost    = total_reduced / nb_pricer

        time_overall = time() - time_beginning
        verbose >= 2 && @printf("Iter %-3d; RMP value: %8.3f;   LG bound: %6.3f,   Pricer: %-8s;   New columns %-2d;   Total columns %-4d;   Aver. red. %5.1f;   Time %5.1f;  Time exact %5.1f; Time RMP: %5.1f  Gap %2.3f\n", 
            iter, rmp_value, lg_bound, "exact", nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, time_pricer, time_rmp, gap)

        if average_reduced_cost > -1. && heuristic_mode
            heuristic_mode=false
            keep_on = true
            verbose >= 2 && println("Switching to exact with validation of routing!")
        end
    end 


    bnp_state.time_rmp      += time_rmp
    bnp_state.time_pricer   += time_pricer

    verbose >= 2 && println("Finished with $lg_bound after $time_overall with $time_pricer for pricer and $time_rmp for RMP")


    try
        if rmp_value > INFEASIBILITY_THRESHOLD
            println("Node $(tree_node.id) infeasible: RMP still uses artificial variables, RMPval: $rmp_value")
            return Inf
        end
        if threshold_ceil(lg_bound) > best_solution-0.01
            return Inf
        end
        branch_or_prune(bnp_state, tree_node, lg_bound, dual_costs)
        return         

    finally
        # unfixing the variables
        unique!(variable_fixed)
        for var in variable_fixed
            unfix(var)
            set_lower_bound(var, 0.0)
            set_upper_bound(var, 1.)
        end

        bnp_state.time_rmp      += time_rmp
        bnp_state.time_pricer   += time_pricer
    end

end






function check_if_all_columns_are_integer(columns)

    for v_subgraph in keys(columns)
        subgraph_is_integer = false
        for column in columns[v_subgraph]
            if value(column.variable) > 0.999
                subgraph_is_integer = true
                break
            end
        end
        if !subgraph_is_integer
            return false
        end
    end

    return true

end






function threshold_ceil(a::Float64; tol::Float64=0.01)
    frac = a - floor(a)
    return frac > tol ? ceil(a) : a
end

