# compact/flow_formulation.jl
# set up and solve the Flow Formulation for VNE.


function solve_flow_formulation(
        instance::Instance;
        time_limit=30.
    )

    time_beginning = time()

    model = Model(CPLEX.Optimizer)
    set_up_formulation!(model, instance)
    set_time_limit_sec(model, time_limit)
    
    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("error! no solution possible...")
        return (
            mapping_cost=Inf,
            solving_time = time()-time_beginning,
            lower_bound = objective_bound(model),
            gap  = relative_gap(model),
            node_count = node_count(model)
        )
    end

    return (
        mapping_cost = objective_value(model),
        solving_time = time()-time_beginning,
        lower_bound = objective_bound(model),
        gap  = relative_gap(model),
        node_count = node_count(model)
    )
end



function solve_flow_formulation_linear(instance::Instance)


    solve_time = 30.

    model = Model(CPLEX.Optimizer)
    set_up_formulation!(model, instance)
    relax_integrality(model)

    set_time_limit_sec(model, solve_time)

    optimize!(model)
    
    println("Resultat : $(objective_value(model))")


end



function set_up_formulation!(model, instance)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    ### Variables

    @variable(model, x[vertices(v_g), vertices(s_g)], binary=true);
    @variable(model, y[edges(v_g), edges(s_dir)], binary=true);

    
    ### Objective

    placement_cost = @expression(model, 
        sum( 
            vn_dem[v_node] * sn_cost[s_node] * x[v_node, s_node] 
                for v_node in vertices(v_g) for s_node in vertices(s_g) 
        )
    )

    routing_cost = @expression(model, sum( ve_dem[src(v_edge), dst(v_edge)] 
                                            * se_cost[src(s_edge), dst(s_edge)] 
                                            * y[v_edge, s_edge]
                                                for v_edge in edges(v_g) 
                                                for s_edge in edges(s_dir) 
    ))

    @objective(model, Min, placement_cost + routing_cost);



    
    ### Pre treatement on capacities 
    
    for v_node in vertices(v_g), s_node in vertices(s_g)
        if vn_dem[v_node] > sn_cap[s_node]
            fix(x[v_node, s_node], 0; force=true)
        end
    end

    for v_edge in edges(v_g), s_edge in edges(s_g)
        if ve_dem[src(v_edge), dst(v_edge)] > se_cap[src(s_edge), dst(s_edge)]
            fix(y[v_edge, s_edge], 0; force=true)
            fix(y[v_edge, Edge(dst(s_edge), src(s_edge))], 0; force=true)
        end
    end

    
    # one substrate node per virtual node
    @constraint(model, [v_node in vertices(v_g)], 
        sum(x[v_node, s_node] for s_node in vertices(s_g)) == 1
    )
  
    # One-to-one node placement
    @constraint(model, [s_node in vertices(s_g)], 
        sum(x[v_node, s_node] for v_node in vertices(v_g)) <= 1
    )


    # Edge capacity constraint (Undirected substrate version)
    @constraint(model, [s_edge in edges(s_g)], 
        sum(
            ve_dem[src(v_edge), dst(v_edge)] * ( y[v_edge, s_edge] + y[v_edge, Edge(dst(s_edge), src(s_edge))]) for v_edge in edges(v_g)) 
            <= se_cap[src(s_edge), dst(s_edge)]
    )
    
    
    # Flow conservation
    @constraint(model, [s_node in vertices(s_g), v_edge in edges(v_g)],
        x[src(v_edge), s_node] - x[dst(v_edge), s_node]
        == sum(y[v_edge, Edge(s_node, s_dst)] for s_dst in outneighbors(s_dir, s_node))
            - sum(y[v_edge, Edge(s_src, s_node)] for s_src in inneighbors(s_dir, s_node))
    )
    
    ## Departure constraints    
    @constraint(model, [s_node in vertices(s_g), v_edge in edges(v_g)],
        sum(y[v_edge, Edge(s_node, s_dst)] for s_dst in outneighbors(s_dir, s_node))
        >= x[src(v_edge), s_node]
    )
    

end



function add_valid_inequalities(model, instance)

    v_g = instance.v_network.graph
    s_g, s_dir = instance.s_network.graph, instance.s_network.directed_graph
    x, y = model[:x], model[:y]

    nb_continuity = 0
    for s_edge_in in edges(s_network_dir)
        for v_edge in edges(v_network)
            s_node = dst(s_edge_in)
            
            if degree(s_network, s_node) < continuity_degree
                if continuity_cap  && (s_network[s_node][:cap]==0)
                    @constraint(model, sum(y[v_edge, s_edge_out] for s_edge_out in get_out_edges(s_network_dir, s_node) ) + x[dst(v_edge), s_node] 
                        >= y[v_edge, s_edge_in] + y[v_edge, Edge(dst(s_edge_in), src(s_edge_in))] )
                    nb_continuity += 1
                elseif !continuity_cap
                    @constraint(model, sum(y[v_edge, s_edge_out] for s_edge_out in get_out_edges(s_network_dir, s_node) ) + x[dst(v_edge), s_node] 
                        >= y[v_edge, s_edge_in] + y[v_edge, Edge(dst(s_edge_in), src(s_edge_in))] )
                    nb_continuity += 1
                end
            end
        end
    end
    
end


function get_solution(model, instance)
    v_g = instance.v_network.graph
    s_g, s_dir = instance.s_network.graph, instance.s_network.directed_graph

end
