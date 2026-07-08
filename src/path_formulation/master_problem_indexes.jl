# src/path_formulation/master_problem.jl
# Master problem and related of the Path Formulation (add reference)

struct Column
    variable
    path
end



function set_up_master_problem!(model::Model, instance::InstanceVNE)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs
    ve_list, se_list = collect(edges(v_g)), collect(edges(s_g))

    # ----- Variables
    @variable(model,  0 <= x[v_node in vertices(v_g), s_node in vertices(s_g)] <= 1)


    # ----- Objective
    placement_cost = @expression(model, sum( vn_dem[v_node] * sn_cost[s_node] * x[v_node, s_node] 
                        for v_node in vertices(v_g), s_node in vertices(s_g) ))
    @objective(model, Min, placement_cost )
    

    # ----- Pretreatment on capacities
    for v_node in vertices(v_g), s_node in vertices(s_g)
        if vn_dem[v_node] > sn_cap[s_node]
            fix(x[v_node, s_node], 0; force=true)
        end
    end


    # ----- Constraints

    # one substrate node per virtual node
    @constraint(model, [v_node in vertices(v_g)], 
        sum(x[v_node, s_node] for s_node in vertices(s_g)) == 1
    )
  
    # One-to-one node placement
    @constraint(model, [s_node in vertices(s_g)], 
        sum(x[v_node, s_node] for v_node in vertices(v_g)) <= 1
    )
    
    # one path per v_edge (relaxed with doi)
    @constraint( model, path_selec[id_ve in 1:ne(v_g)],
        0 >= 1
    )

    # edge capacity
    @constraint( model, capacity_s_edge[id_se in 1:ne(s_g)],
        0 <= se_cap[src(se_list[id_se]), dst(se_list[id_se])]  
    )

    # start of the path
    @constraint( model, start[ve_id in 1:ne(v_g), s_node in vertices(s_g)],
        0 == x[src(ve_list[ve_id]), s_node]
    )
        
    # end of the path
    @constraint( model, terminus[ve_id in 1:ne(v_g), s_node in vertices(s_g)],
        0 == x[dst(ve_list[ve_id]), s_node]
    )

    paths = [ Vector{Column}() for i_edge in 1:ne(v_g)]

    return paths
end




function add_column(model::Model, instance::InstanceVNE, v_edge::Edge, path::Vector{Int}, columns)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    curr_ve_id = instance.v_network.edge_ids[src(v_edge), dst(v_edge)]
    curr_ve_dem = ve_dem[src(v_edge), dst(v_edge)]

    #name_col = "path_$(id_edge)_$(length(columns[i_edge])+1)"
    new_var = @variable(model, lower_bound = 0., upper_bound = 1.0)
    push!(columns[id_edge], Column(new_var, path))

    cost = 0
    for i_node in 1:length(path)-1
        src = path[i_node]
        dst = path[i_node+1]
        cost += curr_ve_dem * se_cost[src, dst]
        se_id = instance.s_network.edge_ids[src, dst]
        set_normalized_coefficient(model[:capacity_s_edge][se_id], new_var, curr_ve_dem)
    end

    set_objective_coefficient(model, new_var, cost)
    set_normalized_coefficient(model[:path_selec][curr_ve_id], new_var, 1)
    
    set_normalized_coefficient(model[:start][ curr_ve_id, path[begin]], new_var, 1)
    set_normalized_coefficient(model[:terminus][ curr_ve_id, path[end]], new_var, 1)  
end


function add_dumb_columns!(model::Model, instance::InstanceVNE)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs
    
    for i_edge in collect(1:ne(v_g))
        i_dumb = 1
        for s_node in vertices(s_g)
            new_var = @variable(model, lower_bound = 0., upper_bound = 1.0)
            set_objective_coefficient(model, new_var, 99999999.)
            set_normalized_coefficient(model[:path_selec][i_edge], new_var, 1)
            set_normalized_coefficient(model[:start][ i_edge, s_node], new_var, 1)
            set_normalized_coefficient(model[:terminus][ i_edge, s_node], new_var, 1)  
            i_dumb+=1
        end
    end
end




struct DualValues
    path_selec::Vector{Float64}
    capacity_edge::Vector{Float64}
    start::Matrix{Float64}
    terminus::Matrix{Float64}
end



# ========= PRICERS PROBLEMS
function solve_pricer(instance, v_edge, dual_costs)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs
    se_ids = instance.s_network.se_ids
    ve_list, se_list, se_dir_list = collect(edges(v_g)), collect(edges(s_g)), collect(edges(s_dir))
    curr_edge_id = ve_list[src(v_edge), dst(v_edge)]

    model_pricer = Model(CPLEX.Optimizer)

    # ----- Variables
    @variable(model_pricer,  x_start[s_node in vertices(s_g)], binary=true)
    @variable(model_pricer,  x_terminus[s_node in vertices(s_g)], binary=true)
    @variable(model_pricer,  y[i_edge in 1:ne(s_dir)], binary=true)

    # ----- Objective
    cost_nodes = @expression(model_pricer, sum( (  - dual_costs.start[curr_edge_id][s_node] ) * x_start[s_node] 
                - dual_costs.terminus[curr_edge_id][s_node] * x_terminus[s_node] for s_node in vertices(s_network)) )
    cost_edges = @expression(model_pricer, sum( ( se_cost[src(se_dir_list[s_edge_id]), dst(se_dir_list[s_edge_id])] - dual_costs.capacity_edge[se_ids[src(se_dir_list[s_edge_id]), dst(se_dir_list[s_edge_id])]] ) * y[s_edge_id] for s_edge_id in 1:ne(s_network_dir)))
    
    @objective(model_pricer, Min, - dual_costs.path_selec[curr_edge_id] + cost_nodes + cost_edges)
    
    
    # ----- Constraints

    # node placement
    @constraint(model, sum( x_start[s_node] for s_node in vertices(s_network)) == 1)
    @constraint(model, sum( x_terminus[s_node] for s_node in vertices(s_network)) == 1)

    # node capacity
    @constraint(model, [s_node in vertices(s_network)],
        x_start[s_node] + x_terminus[s_node] ≤ s_network[s_node][:cap]
    )

    # Flow conservation
    @constraint(model, [s_node in vertices(s_network)],
        sum(  y[s_edge] for s_edge in get_in_edges(s_network_dir, s_node) ) - 
        sum( y[s_edge] for s_edge in get_out_edges(s_network_dir, s_node) ) == 
        x_terminus[s_node] - x_start[s_node] 
    )
    
    # Departure
    @constraint(model, [s_node in vertices(s_network)],
        sum( y[s_edge] for s_edge in get_out_edges(s_network_dir, s_node) ) ≥ x_start[s_node]
    )

    # capacity on nodes
    @constraint(model, [s_node in vertices(s_network)],
        x_start[s_node] + x_terminus[s_node] ≤ 1
    )



    # solve
    set_silent(model)
    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("error! no solution possible...")
        return 
    end

    reduced_cost = objective_value(model)

    if reduced_cost > -0.0001
        return (path_found = nothing, reduced_cost = reduced_cost)
    end

    # Get the solution
    x_start_values = value.(model[:x_start])
    x_terminus_values = value.(model[:x_terminus])
    y_values = value.(model[:y])

    u_start = 0
    u_terminus = 0
    for s_node in  vertices(s_network)
        if x_start_values[s_node] ≥ 0.5
            u_start = s_node
        elseif x_terminus_values[s_node] ≥ 0.5
            u_terminus = s_node
        end
    end

    edges_of_paths = Edge[]
    for  s_edge in edges(s_network_dir)
        if y_values[s_edge] ≥ 0.5
            push!(edges_of_paths, s_edge)
        end
    end

    path = order_path(s_network_dir, edges_of_paths, u_start, u_terminus)

    println("Negative reduced cost! $reduced_cost for $v_edge : $path")


    # return
    return (path_found = path, reduced_cost = reduced_cost)
end




function solve_path_formulation(instance::InstanceVNE)
    model = Model(CPLEX.Optimizer)
    paths = set_up_master_problem!(model, instance);
    add_dumb_columns!(model, instance)
    optimize!(model)
    println("Solved! $(objective_value(model))")
    @time dual(model[:capacity_s_edge][3])
end


