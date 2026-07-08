# path_formulation.jl
# MILP and the polynomial Dijkstra pricers from Moura et al. 2018 for the path formulation



function solve_pricer_exact(instance::Instance, v_edge::Edge, dual_costs::DualValues)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    model_pricer = Model(CPLEX.Optimizer)

    # ----- Variables
    @variable(model_pricer,  x_start[vertices(s_g)], binary=true)
    @variable(model_pricer,  x_terminus[vertices(s_g)], binary=true)
    @variable(model_pricer,  y[edges(s_dir)], binary=true)

    # ----- Objective
    cost_nodes = @expression(model_pricer, sum( - dual_costs.start[v_edge, s_node] * x_start[s_node] 
                - dual_costs.terminus[v_edge, s_node] * x_terminus[s_node] for s_node in vertices(s_g)) )
    cost_edges = @expression(model_pricer, sum( ve_dem[src(v_edge), dst(v_edge)] * (se_cost[src(s_edge), dst(s_edge)] - dual_costs.capacity_edge[s_edge]) * (y[s_edge] + y[Edge(dst(s_edge), src(s_edge))]) for s_edge in edges(s_g) ) )
    
    @objective(model_pricer, Min, - dual_costs.path_selec[v_edge] + cost_nodes + cost_edges)
    
    # ----- Pre-treatment on capacities
    for s_node in vertices(s_g)
        if vn_dem[src(v_edge)] > sn_cap[s_node]
            fix(x_start[s_node], 0.; force=true)
        end
        if vn_dem[dst(v_edge)] > sn_cap[s_node]
            fix(x_terminus[s_node], 0.; force=true)
        end
    end

    for s_edge in edges(s_g)
        if ve_dem[src(v_edge), dst(v_edge)] > se_cap[src(s_edge), dst(s_edge)]
            fix(y[s_edge], 0.; force=true)
            fix(y[Edge(dst(s_edge), src(s_edge))], 0.; force=true)
        end
    end
    
    # ----- Constraints

    # node placement
    @constraint(model_pricer, sum( x_start[s_node] for s_node in vertices(s_g)) == 1)
    @constraint(model_pricer, sum( x_terminus[s_node] for s_node in vertices(s_g)) == 1)

    # one-to-one
    @constraint(model_pricer, [s_node in vertices(s_g)],
        x_start[s_node] + x_terminus[s_node] ≤ 1
    )

    # Flow conservation
    @constraint(model_pricer, [s_node in vertices(s_g)],
        sum( y[Edge(src, s_node)] for src in inneighbors(s_dir, s_node) ) - 
        sum( y[Edge(s_node, dst)] for dst in outneighbors(s_dir, s_node) ) == 
        x_terminus[s_node] - x_start[s_node] 
    )
    
    # Departure
    @constraint(model_pricer, [s_node in vertices(s_g)],
        sum( y[Edge(s_node, dst)] for dst in outneighbors(s_dir, s_node) ) ≥ x_start[s_node]
    )

    # solve
    set_silent(model_pricer)
    optimize!(model_pricer)

    status = primal_status(model_pricer)
    if status != MOI.FEASIBLE_POINT
        println("error! no solution possible...")
        return 
    end

    reduced_cost = objective_value(model_pricer)
    if reduced_cost > -0.0001
        return (path_found = nothing, reduced_cost = reduced_cost)
    end

    # Get the solution
    x_start_values = value.(model_pricer[:x_start])
    x_terminus_values = value.(model_pricer[:x_terminus])
    y_values = value.(model_pricer[:y])

    u_start = 0
    u_terminus = 0
    for s_node in vertices(s_g)
        if x_start_values[s_node] ≥ 0.5
            u_start = s_node
        elseif x_terminus_values[s_node] ≥ 0.5
            u_terminus = s_node
        end
    end

    path = [u_start]
    curr_node = u_start
    while curr_node != u_terminus
        for dst in outneighbors(s_dir, curr_node)
            if y_values[Edge(curr_node, dst)] ≥ 0.5
                push!(path, dst)
                curr_node = dst
                break
            end
        end
    end

    return path, reduced_cost
end


function solve_pricer_dijsktra(instance::Instance, v_edge::Edge, dual_costs::DualValues)
    
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs
    ns = nv(s_g)
    dem_curr_ve = ve_dem[src(v_edge), dst(v_edge)]

    heap = PriorityQueue{Tuple{Int, Int}, Float64}()
    dist = fill(Inf, ns, ns+1)
    prev = zeros(Int, ns, ns+1)
    for s_node_start in vertices(s_g) 
        if vn_dem[src(v_edge)] <= sn_cap[s_node_start]
            heap[(s_node_start, s_node_start)] = - dual_costs.start[v_edge, s_node_start]
            dist[s_node_start, s_node_start] = - dual_costs.start[v_edge, s_node_start]
        end
    end

    while !isempty(heap)

        (u_start, u_curr), curr_dist = dequeue_pair!(heap)
        if u_curr == ns+1
            path = Int[]
            u = ns+1
            while u != u_start
                u = prev[u_start, u]
                push!(path, u)
            end
            reverse!(path)
            reduced_cost = curr_dist - dual_costs.path_selec[v_edge]
            return path, reduced_cost
        end

        for neigh in neighbors(s_g, u_curr)
            if dem_curr_ve <= se_cap[u_curr, neigh]
                curr_cost = 0
                if neigh < u_curr
                    curr_cost = curr_dist + dem_curr_ve * (se_cost[u_curr, neigh] - dual_costs.capacity_edge[Edge(neigh, u_curr)])
                else
                    curr_cost = curr_dist + dem_curr_ve * (se_cost[u_curr, neigh] - dual_costs.capacity_edge[Edge(u_curr, neigh)])
                end
                if curr_cost < dist[u_start, neigh]
                    heap[(u_start, neigh)] = curr_cost
                    prev[u_start, neigh] = u_curr
                    dist[u_start, neigh] = curr_cost
                end
            end
        end

        if u_start != u_curr
            if vn_dem[dst(v_edge)] <= sn_cap[u_curr]
                curr_cost = curr_dist - dual_costs.terminus[v_edge, u_curr]
                if curr_cost < dist[u_start, ns+1]
                    heap[(u_start, ns+1)] = curr_cost
                    prev[u_start, ns+1] = u_curr
                    dist[u_start, ns+1] = curr_cost
                end
            end
        end
    end
    println("No columns found at all!")
    return [], Inf
end

