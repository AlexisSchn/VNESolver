




function solve_branch_price(instance::Instance)

    # 1 ) Set up the decomposition

    time_beginning = time()
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # Compute the partition
    #v_partition = [[1, 2, 3, 11], [5, 6, 9, 11], [4, 7, 12], [8], [10]] # toy instance [[1, 2, 3, 11], [5, 6, 9], [4, 7, 12], [8], [10]] 
    #v_partition = [[15, 16, 17, 18, 19, 24, 26], [1, 2, 20, 21, 22, 27], [9, 10, 11, 12, 13, 15], [4, 5, 6, 14], [3], [7], [8], [23], [25]]
    v_partition = [[i_node] for i_node in 1:nv(v_g)]
    v_partition = [[1, 7, 8, 9, 10, 12, 15], [2, 3, 4, 5, 6], [11, 13, 14]]

    println("Partition: $v_partition")
    v_decomposition = set_up_virtual_decomposition(instance.v_network.graph, v_partition)

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_decomposition.subgraphs)
    println("   and $(length(v_decomposition.cut_edges)) cutting edges")
    println("   and overlapping nodes : $(v_decomposition.overlapping_nodes)")

    
    # RMP
    model_master    = Model(CPLEX.Optimizer)
    set_attribute(model_master, "CPXPARAM_LPMethod", 2)
    set_attribute(model_master, "CPXPARAM_Threads", 1)

    set_silent(model_master)
    print(v_decomposition.nb_appearance_nodes)
    columns         = set_up_master_problem!(model_master, instance, v_decomposition)

    add_dumb_columns!(model_master, v_decomposition)
    add_single_node_columns!(model_master, columns, instance, v_decomposition)

    pricers = Dict{Subgraph, Model}()
    for v_subgraph in v_decomposition.subgraphs
        if length(v_subgraph.nodes) > 1
            pricers[v_subgraph] = set_up_pricer(instance, v_subgraph)
        end
    end


    # 2 ) Solve root node relaxation
    lg_bound, dual_costs = solve_column_generation_root_node2!(model_master, columns, instance, v_decomposition, pricers)
    # Should also get the relaxation itself here! To help with branching. Or I can do it myself here...



    # 4) LOOP TIME
    println("\n\nStarting Branch&Price")

    bnp_state = BranchAndPriceState(instance, v_decomposition, model_master, columns, pricers)
    
    # first branching
    root_node = TreeNode(0, 0, Branching(), 0.) # need to create this
    branch_or_prune(bnp_state, root_node, lg_bound, dual_costs)

    
    time_overall = time() - time_beginning
    best_ub = Inf
    time_max = 500
    nb_nodes_max = 5000
    

    iter = 0
    tree_queue = bnp_state.tree_queue

    while !isempty(tree_queue) && iter < nb_nodes_max && time() - time_beginning < time_max

        current_node = dequeue!(tree_queue)
       
        lg_bound = solve_column_generation_node2!(bnp_state, current_node)

        # prune: collect first, then delete
        to_delete = [node for (node, bound) in tree_queue if threshold_ceil(bound) >= bnp_state.best_solution - 0.01]
        for node in to_delete
            delete!(tree_queue, node)
        end

        best_lb = isempty(tree_queue) ? bnp_state.best_solution : peek(tree_queue)[2]

        gap = (bnp_state.best_solution - best_lb) / best_ub
        time_overall = time() - time_beginning

        @printf("Iter %-3d; UB: %8.1f;  LB: %6.3f,   Nodes left: %3d;   Current node: %3d   Time %5.3f;    Gap %2.3f\n",
            iter, bnp_state.best_solution, best_lb, length(tree_queue), current_node.id, time_overall, gap)

        time_overall = time() - time_beginning
        iter += 1
    end

    println("Finished! Well play")

end