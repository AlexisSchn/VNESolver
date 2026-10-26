# ============================================================================
#   Main function for Branch-and-Price of a VNE instance
#
#   Decomposes the virtual network;
#   Sets the Restricted Master Problem model_master;
#   Sets the Branch&Price
#   Solves root node
#   Then solves the B&B tree
#   
#   @notes this is work-in-progress! 
#   @author Alexis Schneider
# ============================================================================




function solve_branch_price(instance::Instance; params::BranchAndPriceParameters=BranchAndPriceParameters())

    # ===== decomposing the virtual network

    empty!(ROUTING_MODELS)
    time_beginning = time()
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    #v_partition = [[1, 2, 3, 11], [5, 6, 9, 11], [4, 7, 12], [8], [10]] # toy instance [[1, 2, 3, 11], [5, 6, 9], [4, 7, 12], [8], [10]] 
    #v_partition = [[15, 16, 17, 18, 19, 24, 26], [1, 2, 20, 21, 22, 27], [9, 10, 11, 12, 13, 15], [4, 5, 6, 14], [3], [7], [8], [23], [25]]
    v_partition = [[i_node] for i_node in 1:nv(v_g)]
    v_partition = [[1, 7, 8, 9, 10, 12, 15], [2, 3, 4, 5, 6], [11, 13, 14]]
    v_partition = [[8, 9, 10, 12, 15], [2, 3, 4, 5, 6], [1, 6, 7, 11, 13, 14]]

    #v_partition = [[1], [2],[ 3, 4, 5, 6, 7, 8, 9, 25], [9, 10, 11, 12, 15, 16, 17], [13, 14, 15, 18, 19, 20, 21, 22, 23, 24]]
    v_partition = [[7], [3], [19], [18], [1, 2, 25], [4, 5, 6, 8, 9, 25], [10, 11, 12, 15, 16, 17], [13, 14, 15], [15, 20, 21, 22, 23, 24]]
    

    v_decomposition = set_up_virtual_decomposition(instance.v_network.graph, v_partition)
    params.verbose >= 1 && println("Virtual partition: $v_partition")
    params.verbose >= 1 && print_infos_decompo(v_decomposition)


    # ===== builds RMP and B&P utils

    model_master    = Model(CPLEX.Optimizer)
    set_attribute(model_master, "CPXPARAM_LPMethod", 2)
    set_attribute(model_master, "CPXPARAM_Threads", 1)
    set_silent(model_master)
    columns         = set_up_master_problem!(model_master, instance, v_decomposition)
    add_dumb_columns!(model_master, v_decomposition)
    add_single_node_columns!(model_master, columns, instance, v_decomposition)
    pricers         = Dict{Subgraph, Model}()
    for v_subgraph in v_decomposition.subgraphs
        if length(v_subgraph.nodes) > 1
            pricers[v_subgraph] = set_up_pricer(instance, v_subgraph)
        end
    end
    bnp_state       = BranchAndPriceState(instance, v_decomposition, model_master, columns, pricers)





    # ===== solves rootnode

    root_node   = TreeNode(0, 0, Branching(), 0.)
    lg_bound    = solve_column_generation_root_node!(bnp_state)
    optimize!(model_master)
    dual_costs  = DualValues(model_master)
    println("\n\nStarting B&B tree")

    branch_or_prune(bnp_state, root_node, lg_bound, dual_costs)

        
    
    
    # ===== solving B&B

    time_overall    = time() - time_beginning
    best_ub         = Inf
    time_max        = 5000
    nb_nodes_max    = 5000
    iter            = 0
    tree_queue      = bnp_state.tree_queue

    while !isempty(tree_queue) && iter < nb_nodes_max && time() - time_beginning < time_max

        current_node    = dequeue!(tree_queue)
        
        lg_bound        = solve_column_generation_node!(bnp_state, current_node)

        to_delete       = [node for (node, bound) in tree_queue if threshold_ceil(bound) >= bnp_state.best_solution - 0.01]
        for node in to_delete
            delete!(tree_queue, node)
        end
        best_lb         = isempty(tree_queue) ? bnp_state.best_solution : peek(tree_queue)[2]
        gap             = (bnp_state.best_solution - best_lb) / best_ub
        time_overall    = time() - time_beginning
        iter            += 1
        @printf("Iter %-3d; UB: %8.1f;  LB: %6.3f,   Nodes left: %3d;   Current node: %3d   Time %5.1f; Time pricers %5.1f; Time RMP %5.1f   Gap %2.3f\n",
            iter, bnp_state.best_solution, best_lb, length(tree_queue), current_node.id, time_overall, bnp_state.time_pricer, bnp_state.time_rmp, gap)
    end

    println("Finished solving!")

end