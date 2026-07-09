# truc
# You gotta be careful how you mix your dossiers!!

function solve_path_formulation(instance::Instance)
    
    time_beginning = time()
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPXPARAM_LPMethod", 2)
    set_silent(model)
    columns = set_up_master_problem!(model, instance);
    add_dumb_columns!(model, instance)
    add_single_edge_columns!(model, columns, instance)

    time_max        = 1000
    iter_max        = 1000  
    keep_on         = true
    time_overall    = time() - time_beginning
    iter            = 0
    time_pricer     = 0
    time_rmp        = 0

    while keep_on && time_overall < time_max && iter < iter_max
        
        iter += 1        
        keep_on         = false
        total_reduced   = 0
        nb_new_cols     = 0
        
        t=time()
        optimize!(model)
        time_rmp += time()-t

        dual_costs = DualValues(model)
        rmp_value = objective_value(model)

        for v_edge in edges(v_g)
            t=time()
            path, reduced_cost = solve_pricer_dijsktra(instance, v_edge, dual_costs)
            #path, reduced_cost = solve_pricer_exact(instance, v_edge, dual_costs)
            time_pricer += time()-t

            if reduced_cost < - 0.0001
                add_column!(model, columns, instance, v_edge, path)
                keep_on = true
                nb_new_cols += 1
            end
            total_reduced += reduced_cost
        end

        println("Iter $iter; RMP value: $(rmp_value); New columns $nb_new_cols; LG bound: $(rmp_value+total_reduced)")
        time_overall = time() - time_beginning

    end
    
    println("Finished! in $(time()-time_beginning), with $(time_pricer) for pricer and $(time_rmp) for RMP.")

    #=
    # Get the solution for inspection
    println("Node placement:")
    x_values = value.(model[:x])
    for v_node in vertices(v_g)
        println("\tFor vnode $v_node")
        s=0
        for s_node in vertices(s_g)
            if x_values[v_node, s_node] > 0.00001
                s+=x_values[v_node, s_node]
                println("\t\ton $s_node with $(x_values[v_node, s_node])")
                @assert vn_dem[v_node] <= sn_cap[s_node]
            end
        end
        println("\tAt the end, sum $s")
    end
    =#

end



