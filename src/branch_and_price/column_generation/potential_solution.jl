

# After we have found an integer solution, we need to check weither something better can be done 
function cg_after_integer_solution(
        bnp_state::BranchAndPriceState, 
        branching::Branching,
        lg_bound::Float64,
        dual_costs::DualValues
    )

    (; instance, model_master, columns, pricers, v_decomposition, best_solution) = bnp_state

    keep_on = true
    stabilization_coeff = 0.
    rmp_value = objective_value(model_master)


    # need to do another to check that we cant do better
    while keep_on && threshold_ceil(lg_bound) < best_solution - 0.01
        
        keep_on         = false
        total_reduced   = 0
        nb_new_cols     = 0
        old_rmp_value   = rmp_value
        
        optimize!(model_master)
        current_dual_costs = DualValues(model_master)
        stabilize_duals!(dual_costs, current_dual_costs, stabilization_coeff)
        rmp_value = objective_value(model_master)

        if rmp_value < old_rmp_value - 0.01
            if rmp_value < best_solution - 0.01
                return (
                    status = :RMPChanged,
                    lg_bound = lg_bound
                )
            end
        end
        
        for v_subgraph in v_decomposition.subgraphs
            if length(v_subgraph.nodes) == 1
                continue
            end

            submapping, reduced_cost = update_solve_pricer!(pricers[v_subgraph], v_subgraph, dual_costs, instance, branching)

            if reduced_cost < -0.0001 && !isnothing(submapping)
                add_column!(model_master, columns, v_subgraph, submapping, instance)
                keep_on = true
            end

            total_reduced += reduced_cost
        end

        new_lg_bound = rmp_value + total_reduced
        if new_lg_bound > lg_bound
            lg_bound = new_lg_bound
        end

        print(" $lg_bound ")

    end 

    return (
        status = :Valid,
        lg_bound = lg_bound
    )
end



