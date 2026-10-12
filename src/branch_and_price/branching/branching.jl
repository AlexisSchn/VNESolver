



function get_branching(bnp_state::BranchAndPriceState)

    (; instance, model_master, columns, pricers, v_decomposition, best_solution) = bnp_state


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
                    return (
                        status=:NodeBranching,
                        branching_rule=(v_node, best_s_node)
                    )
                end
                continue
            end
        end
    end

    y_values = value.(model_master[:y])

    # I wont care too much about the order for now. 
    # As soon as I find any fractional value, I'll branch
    # off course, a better thing would be to look at 
    for v_edge in v_decomposition.cut_edges
        for s_edge in edges(instance.s_network.directed_graph)
            if y_values[v_edge, s_edge] > 0.01 && y_values[v_edge, s_edge] < 0.99
                return (
                    status=:EdgeBranching,
                    branching_rule=(v_edge, s_edge)
                )
            end
        end
    end


    return (
        status=:IntegerSolution,
        branching_rule=nothing
    )
    
end




function branch_or_prune(
        bnp_state::BranchAndPriceState, 
        tree_node::TreeNode, 
        lg_bound::Float64, 
        dual_costs::DualValues
    )

    (; instance, model_master, columns, pricers, v_decomposition, best_solution, tree_queue) = bnp_state
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    keep_on = true

    while keep_on

        optimize!(model_master)
        rmp_value = objective_value(model_master)

        result_branching = get_branching(bnp_state)

        println(result_branching)

        if result_branching.status == :NodeBranching

            v_node, s_node = result_branching.branching_rule

            # left
            bnp_state.nb_nodes += 1
            branching_left = deepcopy(tree_node.branching)
            branching_left.placement_mandatory[v_node] = s_node
            if v_node ∈ keys(branching_left.placement_allowed)
                delete!(branching_left.placement_allowed, v_node)
            end
            node = TreeNode(bnp_state.nb_nodes, tree_node.depth + 1, branching_left, lg_bound)
            enqueue!(tree_queue, node, node.lower_bound)

            # right
            bnp_state.nb_nodes += 1
            branching_right = deepcopy(tree_node.branching)
            if v_node ∈ keys(branching_right.placement_allowed)
                branching_right.placement_allowed[v_node] = setdiff(branching_right.placement_allowed[v_node], [s_node])
            else
                branching_right.placement_allowed[v_node] = setdiff(vertices(s_g), [s_node])
            end
            node = TreeNode(bnp_state.nb_nodes, tree_node.depth + 1, branching_right, lg_bound)
            enqueue!(tree_queue, node, node.lower_bound)

            return
        
        elseif result_branching.status == :EdgeBranching

            v_edge, s_edge = result_branching.branching_rule
            
            # left
            bnp_state.nb_nodes += 1
            branching_left = deepcopy(tree_node.branching)
            if v_edge ∈ keys(branching_left.routing_cut_edges_forbidden)
                delete!(branching_left.routing_cut_edges_forbidden, v_edge)
            end
            if v_edge ∉ keys(branching_left.routing_cut_edges_mandatory)
                branching_left.routing_cut_edges_mandatory[v_edge] = Edge[]
            end
            push!(branching_left.routing_cut_edges_mandatory[v_edge], s_edge)

            node = TreeNode(nb_node_created, tree_node.depth+1, branching_left, lg_bound)
            enqueue!(tree_queue, node, node.lower_bound)

            # right
            bnp_state.nb_nodes += 1
            branching_right = deepcopy(tree_node.branching)
            if v_edge ∈ keys(branching_right.routing_cut_edges_forbidden)
                push!(branching_right.routing_cut_edges_forbidden[v_edge], s_edge)
            else
                branching_right.routing_cut_edges_forbidden[v_edge] = [s_edge]
            end
            node = TreeNode(nb_node_created, tree_node.depth+1, branching_right, lg_bound)
            enqueue!(tree_queue, node, node.lower_bound)

            return
        
        elseif result_branching.status == :IntegerSolution

            if !check_if_all_columns_are_integer(columns)
                println("For the subgraph, the node placement is integer, but not the routing!")
                error("Edge case not yet tackled")
            end
            
            solution_value = round(Int, rmp_value)

            compare_integer_solution(bnp_state, solution_value)

            result_additional_cg = cg_after_integer_solution(bnp_state, tree_node.branching, lg_bound, dual_costs)

            # CHECK SOME STUFF
            if result_additional_cg.status == :Valid
                return 
            elseif result_additional_cg == :RMPChanged
                lg_bound = result_additional_cg.lg_bound
                println("RMP value has changed: $rmp_value, possible better solution, need to rewind")
                keep_on = true
            end
        end
    end

end
