# src/path_formulation/master_problem.jl
# Master problem and related of the Path Formulation
# Not really the same as Moura 2018

struct Column
    variable::VariableRef
    path::Vector{Int}
end



function set_up_master_problem_false!(model::Model, instance::Instance)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

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

    # one substrate node per virtual node (relaxed)
    @constraint(model, [v_node in vertices(v_g)], 
        sum(x[v_node, s_node] for s_node in vertices(s_g)) >= 1
    )
  
    # One-to-one node placement
    @constraint(model, [s_node in vertices(s_g)], 
        sum(x[v_node, s_node] for v_node in vertices(v_g)) <= 1
    )
    
    # one path per v_edge (relaxed with doi)
    @constraint( model, path_selec[v_edge in edges(v_g)],
        0 >= 1
    )

    # edge capacity 
    @constraint( model, capacity_edge[s_edge in edges(s_g)],
        0 <= se_cap[src(s_edge), dst(s_edge)]  
    )

    # start of the path
    @constraint( model, start[v_edge in edges(v_g), s_node in vertices(s_g)],
        0 <= x[src(v_edge), s_node]
    )
        
    # end of the path
    @constraint( model, terminus[v_edge in edges(v_g), s_node in vertices(s_g)],
        0 <= x[dst(v_edge), s_node]
    )

    columns = Dict(v_edge => Column[] for v_edge in edges(v_g))

    return columns
end


function set_up_master_problem!(model::Model, instance::Instance)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

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

    # one substrate node per virtual node (relaxed)
    @constraint(model, [v_node in vertices(v_g)], 
        sum(x[v_node, s_node] for s_node in vertices(s_g)) == 1
    )
  
    # One-to-one node placement
    @constraint(model, [s_node in vertices(s_g)], 
        sum(x[v_node, s_node] for v_node in vertices(v_g)) <= 1
    )
    
    # one path per v_edge (relaxed with doi)
    @constraint( model, path_selec[v_edge in edges(v_g)],
        0 >= 1
    )

    # edge capacity 
    @constraint( model, capacity_edge[s_edge in edges(s_g)],
        0 <= se_cap[src(s_edge), dst(s_edge)]  
    )

    # start of the path
    @constraint( model, start[v_edge in edges(v_g), s_node in vertices(s_g)],
        0 <= x[src(v_edge), s_node]
    )
        
    # end of the path
    @constraint( model, terminus[v_edge in edges(v_g), s_node in vertices(s_g)],
        0 <= x[dst(v_edge), s_node]
    )

    columns = Dict(v_edge => Column[] for v_edge in edges(v_g))

    return columns
end



function set_up_master_problem_2!(model::Model, instance::Instance)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

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

    # one substrate node per virtual node (relaxed)
    @constraint(model, [v_node in vertices(v_g)], 
        sum(x[v_node, s_node] for s_node in vertices(s_g)) >= 1
    )
  
    # One-to-one node placement
    @constraint(model, [s_node in vertices(s_g)], 
        sum(x[v_node, s_node] for v_node in vertices(v_g)) <= 1
    )
    
    # one path per v_edge (relaxed with doi)
    @constraint( model, path_selec[v_edge in edges(v_g)],
        0 >= 1
    )

    # edge capacity 
    @constraint( model, capacity_edge[s_edge in edges(s_g)],
        0 <= se_cap[src(s_edge), dst(s_edge)]  
    )

    # start of the path
    @constraint( model, start[v_edge in edges(v_g), s_node in vertices(s_g)],
        0 >= x[src(v_edge), s_node]
    )
        
    # end of the path
    @constraint( model, terminus[v_edge in edges(v_g), s_node in vertices(s_g)],
        0 >= x[dst(v_edge), s_node]
    )

    columns = Dict(v_edge => Column[] for v_edge in edges(v_g))

    return columns
end


function add_column!(model::Model, columns, instance::Instance, v_edge::Edge, path::Vector{Int})

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    curr_ve_dem = ve_dem[src(v_edge), dst(v_edge)]

    new_var = @variable(model, lower_bound = 0., upper_bound = 1.0)
    push!(columns[v_edge], Column(new_var, path))

    cost = 0
    for i_node in 1:length(path)-1
        src = path[i_node]
        dst = path[i_node+1]
        if dst < src
            src, dst = dst, src
        end
        cost += curr_ve_dem * se_cost[src, dst]
        set_normalized_coefficient(model[:capacity_edge][Edge(src, dst)], new_var, curr_ve_dem)
    end

    set_objective_coefficient(model, new_var, cost)
    set_normalized_coefficient(model[:path_selec][v_edge], new_var, 1)
    
    set_normalized_coefficient(model[:start][ v_edge, path[begin]], new_var, 1)
    set_normalized_coefficient(model[:terminus][ v_edge, path[end]], new_var, 1)  
end


function add_dumb_columns!(model::Model, instance::Instance)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs
    
    for v_edge in edges(v_g)
        new_var = @variable(model, lower_bound = 0., upper_bound = 1.0)
        set_objective_coefficient(model, new_var, 10e5)
        set_normalized_coefficient(model[:path_selec][v_edge], new_var, 1)
    end
end


function add_single_edge_columns!(model::Model, columns, instance::Instance)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    nb_cols_init = 0
    for v_edge in edges(v_g)
        for s_edge in edges(s_dir)
            if se_cap[src(s_edge), dst(s_edge)] >= ve_dem[src(v_edge), dst(v_edge)]
                if sn_cap[src(s_edge)] >= vn_dem[src(v_edge)] && sn_cap[dst(s_edge)] >= vn_dem[dst(v_edge)]
                    add_column!(model, columns, instance, v_edge, [src(s_edge), dst(s_edge)])
                    nb_cols_init += 1
                end
            end
        end
    end
    println("Found $nb_cols_init single-edge paths for initialization")
end


struct DualValues
    path_selec::JuMP.Containers.DenseAxisArray{Float64, 1}
    capacity_edge::JuMP.Containers.DenseAxisArray{Float64, 1}
    start::JuMP.Containers.DenseAxisArray{Float64, 2}
    terminus::JuMP.Containers.DenseAxisArray{Float64, 2}
end

function DualValues(model::Model)
    return DualValues(
        dual.(model[:path_selec]),
        dual.(model[:capacity_edge]),
        dual.(model[:start]),
        dual.(model[:terminus])
    )
end



