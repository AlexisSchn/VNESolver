

function solve_column_generation_node!(
        model_master::Model, 
        columns, 
        instance::Instance, 
        v_decomposition::OverlappingVirtualDecomposition, 
        pricers, 
        tree_node::TreeNode,
        tree_queue::PriorityQueue{TreeNode, Float64},
        nb_node_created::Int,
        best_solution
    )

    time_beginning = time()
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs
    

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
    for v_edge in keys(tree_node.branching.routing_mandatory)
        s_edge = tree_node.branching.routing_mandatory[v_edge]
        var = model_master[:y][v_edge, s_edge]
        fix(var, 1.; force=true)
        push!(variable_fixed, var)
    end
    for v_edge in keys(tree_node.branching.routing_forbidden)
        for s_edge in tree_node.branching.routing_forbidden[v_edge]
            var = model_master[:y][v_edge, s_edge]
            fix(var, 0.; force=true)
            push!(variable_fixed, var)
        end
    end


    # =============== COLUMN GENERATION

    optimize!(model_master)

    dual_costs = DualValues(model_master)


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



    # greedy pricers
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

        modified_se_cost = zeros(Float64, nv(s_g), nv(s_g))
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
        
            submapping, reduced_cost = solve_greedy_pricer_branching(instance, v_subgraph, current_dual_costs, shortest_paths.dists, modified_se_cost, tree_node.branching, nb_greedy = 50, time_max = 10)
                
        
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


    # exact pricers
    keep_on=true
    optimize!(model_master)
    dual_costs = DualValues(model_master)
    while keep_on && time_overall < time_max && iter < nb_iter_max && gap > gap_min && threshold_ceil(lg_bound) < best_solution - 0.01

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


    # POST - CG : solution or branching decision

    optimize!(model_master)
    rmp_value = objective_value(model_master)

    if objective_value(model_master) > 10e5 # why this value here? Idk what you are trying to achieve...
        println("Large value for master problem $(objective_value(model_master)), but got $lg_bound LBG")
        println("THIS IS NOT TAKEN INTO ACCOUNT YET!!")
        error("Sorry we gotta stop")
    end

    try
        if lg_bound > 10e5
            println("Node unfeasible")
            return (lg_bound, Inf)
        end

        if threshold_ceil(lg_bound) > best_solution-0.01
            return (lg_bound, Inf)
        end


        columns_all_integer = check_if_all_columns_are_integer(columns)
        if !columns_all_integer

            (v_node, s_node) = get_x_branching(model_master, instance, v_decomposition, columns)

            # left
            branching_left = deepcopy(tree_node.branching)
            branching_left.placement_mandatory[v_node] = s_node
            if v_node ∈ keys(branching_left.placement_allowed)
                delete!(branching_left.placement_allowed, v_node)
            end
            nb_node_created+=1
            node = TreeNode(nb_node_created, tree_node.depth+1, branching_left, lg_bound)
            enqueue!(tree_queue, node, node.lower_bound)

            # right
            branching_right = deepcopy(tree_node.branching)
            if v_node ∈ keys(branching_right.placement_allowed)
                branching_right.placement_allowed[v_node] = setdiff(branching_right.placement_allowed[v_node], [s_node])
            else
                branching_right.placement_allowed[v_node] = setdiff(vertices(s_g), [s_node])
            end
            nb_node_created+=1
            node = TreeNode(nb_node_created, tree_node.depth+1, branching_right, lg_bound)
            enqueue!(tree_queue, node, node.lower_bound)

        else

            (v_edge, s_edge) = get_y_branching(model_master, instance, v_decomposition)
            
            if isnothing(v_edge) # Integer solution!
                println("New solution found! value : $rmp_value and lgb : $lg_bound")
                
                solution_value = round(Int, rmp_value)
                best_sol = minimum([best_solution, solution_value])
                
                # need to do another to check that we cant do better
                while keep_on && threshold_ceil(lg_bound) < best_solution - 0.01
                    
                    old_rmp_value = rmp_value
                    keep_on         = false
                    total_reduced   = 0
                    nb_new_cols     = 0
                    
                    t=time()

                    optimize!(model_master)
                    time_rmp += time()-t
                    current_dual_costs = DualValues(model_master)
                    stabilize_duals!(dual_costs, current_dual_costs, stabilization_coeff)
                    rmp_value = objective_value(model_master)

                    if rmp_value < rmp_value - 0.01
                        println("OMG! THE VALUE HAS CHANGED!!")
                    end

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

                    print(" $lg_bound ")
                    #@printf("Iter %-3d; RMP value: %8.3f;   LG bound: %6.3f,   Pricer: %-8s;   New columns %-2d;   Total columns %-4d;   Aver. red. %5.3f;   Time %5.3f;    Gap %2.3f\n", 
                    #    iter, rmp_value, lg_bound, "milp", nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, gap)
                end 

                println("AT THE END : $lg_bound")

                return (lg_bound, solution_value)

            else # need to branch

                # left
                branching_left = deepcopy(tree_node.branching)
                branching_left.routing_mandatory[v_edge] = s_edge
                if v_edge ∈ keys(branching_left.routing_forbidden)
                    delete!(branching_left.routing_forbidden, v_edge)
                end
                nb_node_created+=1
                node = TreeNode(nb_node_created, tree_node.depth+1, branching_left, lg_bound)
                enqueue!(tree_queue, node, node.lower_bound)

                # right
                branching_right = deepcopy(tree_node.branching)
                if v_edge ∈ keys(branching_right.routing_forbidden)
                    push!(branching_right.routing_forbidden[v_edge], s_edge)
                else
                    branching_right.routing_forbidden[v_edge] = [s_edge]
                end
                nb_node_created+=1
                node = TreeNode(nb_node_created, tree_node.depth+1, branching_right, lg_bound)
                enqueue!(tree_queue, node, node.lower_bound)
            end
        end

        return (lg_bound, Inf)

    finally

        # Unfixing the stuff
        unique!(variable_fixed)
        for var in variable_fixed
            unfix(var)
            set_lower_bound(var, 0.0)
            set_upper_bound(var, 1.)
        end

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


function get_x_branching(model_master::Model, instance::Instance, v_decomposition::OverlappingVirtualDecomposition, columns)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs


    # First, get the nodes order
    overlapping_in_vg = Int[]
    non_overlapping_in_vg = Int[]
    for v in vertices(v_g)
        if v ∈ v_decomposition.overlapping_nodes
            push!(overlapping_in_vg, v)
        else
            push!(non_overlapping_in_vg, v)
        end
    end
    sort!(overlapping_in_vg, by = v -> degree(v_g, v), rev = true)
    sort!(non_overlapping_in_vg, by = v -> degree(v_g, v), rev = true)
    branching_order = vcat(overlapping_in_vg, non_overlapping_in_vg)
    
    #println("The branching order is: $branching_order")


    # now, look at the solution 
    for v_node in branching_order
        for v_subgraph in v_decomposition.subgraphs 
            if v_node ∈ v_subgraph.nodes
                # Let's look at the fractional values
                placement_s_nodes = zeros(nv(s_g))
                for column in columns[v_subgraph]
                    val = value(column.variable)
                    if val > 0.00001
                        placement_s_nodes[column.submapping.node_placement[v_node]] += val
                    end
                end
                #println(placement_s_nodes)
                if sum(placement_s_nodes) < 0.99
                    println("Objective rn: $(objective_value(model_master))")
                    error("AHh THE SUM IS NOT 1")
                end
                if maximum(placement_s_nodes) < 0.99
                    best_s_node = argmax(placement_s_nodes)
                    #println("We will branch $v_node on $best_s_node, which is at $(placement_s_nodes[best_s_node])")
                    return (v_node, best_s_node)
                end
                continue
            end
        end
    end

    println("WHAT THE FUUUUCK")
    println("The columns are fractional but somehow all the placement is integer?")
    return (0, 0)
end


function get_y_branching(model_master::Model, instance::Instance, v_decomposition::OverlappingVirtualDecomposition)

    y_values = value.(model_master[:y])

    # I wont care too much about the order for now. 
    # As soon as I find any fractional value, I'll branch
    # off course, a better thing would be to lookl 
    for v_edge in v_decomposition.cut_edges
        for s_edge in edges(instance.s_network.directed_graph)
            if y_values[v_edge, s_edge] > 0.01 && y_values[v_edge, s_edge] < 0.99
                return (v_edge, s_edge)
            end
        end
    end

    return (nothing, nothing)
end



function threshold_ceil(a::Float64; tol::Float64=0.01)
    frac = a - floor(a)
    return frac > tol ? ceil(a) : a
end

