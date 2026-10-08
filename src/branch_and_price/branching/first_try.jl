



function first_idea(model_master::Model, instance::Instance, v_decomposition::OverlappingVirtualDecomposition, columns)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs


    # First, get the nodes order
    overlapping_in_vg = Int[]
    non_overlapping_in_vg = Int[]
    for v in vertices(v_g)
        if v ∈ v_decomposition.overlapping_nodes
            push!(overlapping_in_vg, v)
        else
            push!(non_overlapping_in_vg, v)
        end
    end
    sort!(overlapping_in_vg, by = v -> degree(v_g, v), rev = true)
    sort!(non_overlapping_in_vg, by = v -> degree(v_g, v), rev = true)
    branching_order = vcat(overlapping_in_vg, non_overlapping_in_vg)
    
    #println("The branching order is: $branching_order")


    # now, look at the solution 
    for v_node in branching_order
        for v_subgraph in v_decomposition.subgraphs 
            if v_node ∈ v_subgraph.nodes
                # Let's look at the fractional values
                placement_s_nodes = zeros(nv(s_g))
                for column in columns[v_subgraph]
                    val = value(column.variable)
                    if val > 0.00001
                        placement_s_nodes[column.submapping.node_placement[v_node]] += val
                    end
                end
                #println(placement_s_nodes)
                if sum(placement_s_nodes) < 0.99
                    println("Objective rn: $(objective_value(model_master))")
                    error("AHh THE SUM IS NOT 1")
                end
                if maximum(placement_s_nodes) < 0.99
                    best_s_node = argmax(placement_s_nodes)
                    #println("We will branch $v_node on $best_s_node, which is at $(placement_s_nodes[best_s_node])")
                    return (v_node, best_s_node)
                end
                continue
            end
        end
    end

    return (0, 0)
end


function second_try(model_master::Model, instance::Instance, v_decomposition::OverlappingVirtualDecomposition, columns, rmp_value)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs


    # First, get the nodes order
    overlapping_in_vg = Int[]
    non_overlapping_in_vg = Int[]
    for v in vertices(v_g)
        if v ∈ v_decomposition.overlapping_nodes
            push!(overlapping_in_vg, v)
        else
            push!(non_overlapping_in_vg, v)
        end
    end
    sort!(overlapping_in_vg, by = v -> degree(v_g, v), rev = true)
    sort!(non_overlapping_in_vg, by = v -> degree(v_g, v), rev = true)
    branching_order = vcat(overlapping_in_vg, non_overlapping_in_vg)
    
    #println("The branching order is: $branching_order")


    # now, look at the solution 
    for v_node in branching_order
        for v_subgraph in v_decomposition.subgraphs 
            if v_node ∈ v_subgraph.nodes
                # Let's look at the fractional values
                placement_s_nodes = zeros(nv(s_g))
                for column in columns[v_subgraph]
                    val = value(column.variable)
                    if val > 0.00001
                        placement_s_nodes[column.submapping.node_placement[v_node]] += val
                    end
                end
                #println(placement_s_nodes)
                if sum(placement_s_nodes) < 0.999
                    println("Objective rn: $(objective_value(model_master))")
                    error("AHh THE SUM IS NOT 1")
                end
                if maximum(placement_s_nodes) < 0.999 # otherwise, its integer!
                    best_s_node = argmax(placement_s_nodes)
                    #println("We will branch $v_node on $best_s_node, which is at $(placement_s_nodes[best_s_node])")
                    return (v_node, best_s_node)
                end
                continue
            end
        end
    end

    # if we are here, it means that all 
    println("Well!!! It seems that everything is integer... Congrats...")
    return (0, 0)
end