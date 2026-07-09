


function solve_subgraph_decomposition(instance::Instance; parameters::SubgraphDecompositionParameters=SubgraphDecompositionParameters())

    time_beginning = time()
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    # Compute the partition
    #v_partition = [[i] for i in 1:nv(v_g)]
    v_partition = star_partition(v_g)
    println("Partition: $v_partition")

    # Setting the decomposition
    v_subgraphs = []
    cut_edges = Edge[]
    for v_nodes in v_partition
        sort!(v_nodes)
        idx_of_nodes = zeros(Int, nv(v_g))
        for (i_node, v_node) in enumerate(v_nodes) 
            idx_of_nodes[v_node] = i_node
        end

        subgraph_edges = Edge[]
        cut_edges_with_src = [Edge[] for _ in v_nodes]
        cut_edges_with_dst = [Edge[] for _ in v_nodes]

        for v_edge in edges(v_g)
            if src(v_edge) ∈ v_nodes
                if dst(v_edge) ∈ v_nodes
                    push!(subgraph_edges, v_edge)
                else
                    push!(cut_edges, v_edge)
                    v_source = findfirst(==(src(v_edge)), v_nodes) 
                    push!(cut_edges_with_src[v_source], v_edge)  
                end
            end
            if dst(v_edge) ∈ v_nodes && src(v_edge) ∉ v_nodes
                v_dst = findfirst(==(dst(v_edge)), v_nodes) 
                push!(cut_edges_with_dst[v_dst], v_edge)  
            end
        end
        push!(v_subgraphs, Subgraph(v_nodes, subgraph_edges, cut_edges_with_src, cut_edges_with_dst, idx_of_nodes))
    end
    v_decomposition = VirtualDecomposition(v_subgraphs, cut_edges)
    
    # RMP
    model_master    = Model(CPLEX.Optimizer)
    set_attribute(model_master, "CPXPARAM_LPMethod", 2)
    set_silent(model_master)
    columns         = set_up_master_problem!(model_master, instance, v_decomposition)

    add_dumb_columns!(model_master, v_decomposition)
    add_single_node_columns!(model_master, columns, instance, v_decomposition)
    add_greedy_columns!(model_master, columns, instance, v_decomposition)

    # Column generation
    pricers = Dict{Subgraph, Model}()
    for v_subgraph in v_subgraphs
        pricers[v_subgraph] = set_up_pricer(instance, v_subgraph)
    end
    time_max            = parameters.time_max
    nb_iter_max         = parameters.nb_iter_max
    gap_min             = parameters.gap_min
    keep_on             = true

    time_overall        = time() - time_beginning
    iter                = 0
    gap                 = 1.
    time_pricer         = 0
    time_rmp            = 0
    rmp_value, lg_bound = 0.,0.
    nb_total_columns    = 0
    nb_pricer           = 0
    for v_subgraph in v_subgraphs
        if length(v_subgraph.nodes) > 1
            nb_pricer += 1
        end
    end
    pricer              = "greedy"

    while keep_on && time_overall < time_max && iter < nb_iter_max && gap > gap_min
        
        iter += 1        
        keep_on         = false
        total_reduced   = 0
        nb_new_cols     = 0
        
        t=time()

        optimize!(model_master)
        time_rmp += time()-t

        dual_costs = DualValues(model_master)
        rmp_value = objective_value(model_master)

        for v_subgraph in v_subgraphs
            if length(v_subgraph.nodes) == 1
                continue
            end

            t=time()
            
            if pricer=="greedy"
                submapping, reduced_cost = solve_greedy_pricer(instance, v_subgraph, dual_costs, nb_greedy = 100, time_max = 10)
            elseif pricer=="exact"
                submapping, reduced_cost = update_solve_pricer!(pricers[v_subgraph], v_subgraph, dual_costs, instance)
            end
            time_pricer += time()-t

            if reduced_cost < -0.0001 && !isnothing(submapping)
                add_column!(model_master, columns, v_subgraph, submapping, instance)
                keep_on = true
                nb_new_cols += 1
                nb_total_columns += 1
            end

            total_reduced += reduced_cost
        end

        if pricer=="exact"
            new_lg_bound = rmp_value + total_reduced
            if new_lg_bound > lg_bound
                lg_bound = new_lg_bound
            end
        end
        gap = (rmp_value-lg_bound)/rmp_value
        average_reduced_cost = total_reduced / nb_pricer
        if average_reduced_cost > -3. && pricer=="greedy" && rmp_value < 10000
            pricer="exact"
            println("Switching to exact pricers!")
        end

        time_overall = time() - time_beginning
        @printf("Iter %-3d; RMP value: %8.3f;   LG bound: %6.3f,   Pricer: %-8s;   New columns %-2d;   Total columns %-4d;   Aver. red. %5.3f;   Time %5.3f;    Gap %2.3f\n", 
            iter, rmp_value, lg_bound, pricer, nb_new_cols, nb_total_columns, average_reduced_cost, time_overall, gap)

    end
    
    result = SubgraphDecompositionResult(
        instance.v_network.name,
        instance.s_network.name,
        rmp_value, 
        lg_bound,
        gap,
        nb_total_columns,
        iter,
        time()-time_beginning
    )
    println("Result: $result")
    return result
end




function add_greedy_columns!(model_master::Model, columns, instance::Instance, v_decomposition::VirtualDecomposition)
    empty_duals = zero_duals(model_master)
    for v_subgraph in v_decomposition.subgraphs
        if length(v_subgraph.nodes) <= 2
            continue
        end
        for i in 1:50
            submapping, reduced_cost = solve_greedy_pricer(instance, v_subgraph, empty_duals, nb_greedy = 3, time_max = 10)
            if isnothing(submapping)
                continue
            end
            #println("Found a nice mapping $reduced_cost")
            #=
            if column_already_there(columns, submapping, v_subgraph)
                println("already there!")
                continue
            end
            =#
            add_column!(model_master, columns, v_subgraph, submapping, instance)
        end
    end
end




function star_partition(graph::Graph)

    copy_v_network = copy(graph)
    keep_on = true

    real_indices = collect(1:nv(graph)) 
    # Graphs.jl is completly stupid when it comes to removing a node. It swaps it with the last node and then removes it. be careful... 

    nodes_in_no_subgraphs = collect(1:nv(graph))
    v_node_partitionning = Vector{Vector{Int}}()
    possible_centers = collect(1:nv(graph))

    while keep_on

        # Get node with max degree
        nodes_degrees = [degree(copy_v_network, v_node) for v_node in vertices(copy_v_network)]
        node_sorted_degree = sortperm(nodes_degrees, rev=true)
        center = 0
        for node in node_sorted_degree
            if real_indices[node] ∈ possible_centers
                center = node
                break
            end
        end
        if degree(copy_v_network, center) <= 1
            break
        end

        new_part = [ real_indices[center] ]
        for neigh in neighbors(copy_v_network, center)
            #println("Neighbor: $neigh, aka $(real_indices[neigh])")
            push!(new_part, real_indices[neigh])
        end

        for v_node in new_part
            filter!(x -> x != v_node, possible_centers)
        end
    
        sort!(new_part)
        push!(v_node_partitionning, new_part)

        for v_node in new_part
            if v_node ∈ nodes_in_no_subgraphs
                filter!(x -> x != v_node, nodes_in_no_subgraphs)
            end
        end

        while length(neighbors(copy_v_network, center)) > 0
            neigh = neighbors(copy_v_network, center)[1]
            rem_vertex!(copy_v_network, neigh)

            real_indices[neigh] = real_indices[length(real_indices)]
            if center == length(real_indices)
                center = neigh
            end
            deleteat!(real_indices, length(real_indices))


        end

        rem_vertex!(copy_v_network, center)
        real_indices[center] = real_indices[length(real_indices)]
        deleteat!(real_indices, length(real_indices))

    end


    for v_node_left in nodes_in_no_subgraphs
        push!(v_node_partitionning, [v_node_left])
    end

    return v_node_partitionning

end

