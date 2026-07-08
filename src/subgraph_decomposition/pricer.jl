# Pricers...



function set_up_pricer(instance::Instance, v_subgraph::Subgraph)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    model = Model(CPLEX.Optimizer)

    
    ### Variables
    @variable(model, x[v_subgraph.nodes, vertices(s_g)], binary=true);
    @variable(model, y[v_subgraph.edges, edges(s_dir)], binary=true);
    
    ### Objective: empty for nows
    @objective(model, Min, 0.);

    
    ### Pre treatement on capacities 
    
    for v_node in v_subgraph.nodes, s_node in vertices(s_g)
        if vn_dem[v_node] > sn_cap[s_node]
            fix(x[v_node, s_node], 0; force=true)
        end
    end

    for v_edge in v_subgraph.edges, s_edge in edges(s_g)
        if ve_dem[src(v_edge), dst(v_edge)] > se_cap[src(s_edge), dst(s_edge)]
            fix(y[v_edge, s_edge], 0; force=true)
            fix(y[v_edge, Edge(dst(s_edge), src(s_edge))], 0; force=true)
        end
    end

    
    # one substrate node per virtual node
    @constraint(model, [v_node in v_subgraph.nodes], 
        sum(x[v_node, s_node] for s_node in vertices(s_g)) == 1
    )
  
    # One-to-one node placement
    @constraint(model, [s_node in vertices(s_g)], 
        sum(x[v_node, s_node] for v_node in v_subgraph.nodes) <= 1
    )


    # Edge capacity constraint (Undirected substrate version)
    @constraint(model, [s_edge in edges(s_g)], 
        sum(
            ve_dem[src(v_edge), dst(v_edge)] * ( y[v_edge, s_edge] + y[v_edge, Edge(dst(s_edge), src(s_edge))]) for v_edge in v_subgraph.edges) 
            <= se_cap[src(s_edge), dst(s_edge)]
    )
    
    
    # Flow conservation
    @constraint(model, [s_node in vertices(s_g), v_edge in v_subgraph.edges],
        x[src(v_edge), s_node] - x[dst(v_edge), s_node]
        == sum(y[v_edge, Edge(s_node, s_dst)] for s_dst in outneighbors(s_dir, s_node))
            - sum(y[v_edge, Edge(s_src, s_node)] for s_src in inneighbors(s_dir, s_node))
    )
    
    ## Departure constraints    
    @constraint(model, [s_node in vertices(s_g), v_edge in v_subgraph.edges],
        sum(y[v_edge, Edge(s_node, s_dst)] for s_dst in outneighbors(s_dir, s_node))
        >= x[src(v_edge), s_node]
    )
    
    return model
end


function update_solve_pricer!(model::Model, v_subgraph::Subgraph, duals::DualValues, instance::Instance)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # updating costs
    placement_cost = @expression(model, 
        sum( ( sn_cost[s_node] * vn_dem[v_node] - duals.node_1t1[s_node] ) * model[:x][v_node, s_node] 
            for v_node in v_subgraph.nodes for s_node in vertices(s_g) ))
  
               
    routing_cost = @expression(model, sum( 
        ( se_cost[src(s_edge), dst(s_edge)] - duals.edge_capacity[s_edge] ) * ve_dem[src(v_edge), dst(v_edge)]
        *  (model[:y][v_edge, s_edge] + model[:y][v_edge, Edge(dst(s_edge), src(s_edge))])
                for v_edge in v_subgraph.edges for s_edge in edges(s_g) ))


            
    # flow conservation
    flow_conservation_cost = @expression( model, 
        sum(
            - sum(duals.flow_conservation[s_node, cut_edge] * model[:x][src(cut_edge), s_node] for cut_edge in v_subgraph.cut_edges_with_src[i_node])
            + sum(duals.flow_conservation[s_node, cut_edge] * model[:x][dst(cut_edge), s_node]  for cut_edge in v_subgraph.cut_edges_with_dst[i_node])
                for s_node in vertices(s_g), i_node in 1:length(v_subgraph.nodes) )
    )
    
    # departure
    departure_costs = @expression( model,
        - sum( duals.flow_departure[s_node, cut_edge] * model[:x][src(cut_edge), s_node] 
            for s_node in vertices(s_g), i_node in 1:length(v_subgraph.nodes), cut_edge in v_subgraph.cut_edges_with_src[i_node])
    )


    
    @objective(model, Min, 
        - duals.submapping_selection[v_subgraph]
        + placement_cost + routing_cost 
        + flow_conservation_cost 
        + departure_costs
    );

    # solving
    set_silent(model)
    optimize!(model)

    status = primal_status(model)
    if (status != MOI.FEASIBLE_POINT)
        println("Infeasible subproblem... $status")
        return nothing, Inf
    end

    reduced_cost = objective_value(model) 
    if reduced_cost > -0.000001
        return nothing, 0
    end


    # returning the solution
    x_values = value.(model[:x])
    y_values = value.(model[:y])

    node_placement = zeros(Int, length(v_subgraph.nodes))
    for (i_node, v_node) in enumerate(v_subgraph.nodes)
        for s_node in vertices(s_g)
            if x_values[v_node, s_node] > 0.1
                node_placement[i_node] = s_node
            end
        end
    end

    edge_routing = Vector{Vector{Int}}()
    for v_edge in v_subgraph.edges
        idx_src, idx_dst = findfirst(==(src(v_edge)), v_subgraph.nodes), findfirst(==(dst(v_edge)), v_subgraph.nodes) 
        u_start, u_terminus = node_placement[idx_src], node_placement[idx_dst]
        path = [u_start]
        curr_node = u_start
        while curr_node != u_terminus
            for dst in outneighbors(s_dir, curr_node)
                if y_values[v_edge, Edge(curr_node, dst)] ≥ 0.5 
                    push!(path, dst)
                    curr_node = dst                
                    break
                end
            end
        end
        push!(edge_routing, path)
    end

    return Mapping(node_placement, edge_routing), reduced_cost

end
