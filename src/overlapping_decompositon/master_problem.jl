

struct Column
    variable
    submapping::Mapping
end


function set_up_master_problem!(model_master::Model, instance::Instance, v_decomposition::OverlappingVirtualDecomposition)

    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs
    
        
    # ----- Variables
    @variable(model_master, 0. <= y[v_decomposition.cut_edges, edges(s_dir)] <= 1.);

    
    # ----- Variables
    master_routing_cost = @expression(model_master, sum( ve_dem[src(v_edge), dst(v_edge)] * se_cost[src(s_edge), dst(s_edge)] * y[v_edge, s_edge] 
                        for v_edge in v_decomposition.cut_edges, s_edge in edges(s_dir) ))
    @objective(model_master, Min, master_routing_cost )


    # ----- Pre-treatment on capacities
    for v_edge in v_decomposition.cut_edges, s_edge in edges(s_g)
        if ve_dem[src(v_edge), dst(v_edge)] > se_cap[src(s_edge), dst(s_edge)]
            fix(y[v_edge, s_edge], 0; force=true)
            fix(y[v_edge, Edge(dst(s_edge), src(s_edge))], 0; force=true)
        end
    end

    # one substrate submapping per virtual subgraph
    @constraint(model_master, submapping_selection[v_subgraph in v_decomposition.subgraphs], 
        0 >= 1
    )

    # overlapping nodes
    @constraint(model_master, overlapping[s_node in vertices(s_g), v_node in v_decomposition.overlapping_nodes, i_subgraph in 1:(v_decomposition.nb_appearance_nodes[v_node])],
        0 <= 0
    )
  
    # One-to-one node placement
    @constraint(model_master, node_1t1[s_node in vertices(s_g)], 
        0 <= 1
    )
    
    # Edge capacity constraint (Undirected substrate version)
    @constraint(model_master, edge_capacity[s_edge in edges(s_g)], 
        sum(
            0 + ve_dem[src(v_edge), dst(v_edge)] * ( y[v_edge, s_edge] + y[v_edge, Edge(dst(s_edge), src(s_edge))]) for v_edge in v_decomposition.cut_edges) 
            <= se_cap[src(s_edge), dst(s_edge)]
    )
    
    # Flow conservation
    @constraint(model_master, flow_conservation[s_node in vertices(s_g), v_edge in v_decomposition.cut_edges],
        0
        == sum(y[v_edge, Edge(s_node, s_dst)] for s_dst in outneighbors(s_dir, s_node))
            - sum(y[v_edge, Edge(s_src, s_node)] for s_src in inneighbors(s_dir, s_node))
    )
    
    ## Departure constraints    
    @constraint(model_master, flow_departure[s_node in vertices(s_g), v_edge in v_decomposition.cut_edges],
        0 <= sum(y[v_edge, Edge(s_node, s_dst)] for s_dst in outneighbors(s_dir, s_node))
    )
    

    columns = Dict(v_subgraph => Column[] for v_subgraph in v_decomposition.subgraphs)

    return columns
end



function add_column!(model_master::Model, columns, v_subgraph::Subgraph, submapping::Mapping, instance::Instance)
      
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs
    
    new_var = @variable(model_master, lower_bound = 0., upper_bound = 1.0)
    push!(columns[v_subgraph], Column(new_var, submapping))

    cost = 0
    for v_node in v_subgraph.nodes 
        s_node = submapping.node_placement[v_node]

        cost += vn_dem[v_node] / v_subgraph.nb_appearance_nodes[v_node] * sn_cost[s_node] #! Very important, be careful here! This changes everything!

        set_normalized_coefficient(model_master[:node_1t1][s_node], new_var, 1/v_subgraph.nb_appearance_nodes[v_node])
        for cut_edge in v_subgraph.cut_edges_src[v_node]
            set_normalized_coefficient(model_master[:flow_conservation][s_node, cut_edge], new_var, 1/v_subgraph.nb_appearance_nodes[v_node] )
            set_normalized_coefficient(model_master[:flow_departure][s_node, cut_edge], new_var, 1/v_subgraph.nb_appearance_nodes[v_node] )
        end
        for cut_edge in v_subgraph.cut_edges_dst[v_node]
            set_normalized_coefficient(model_master[:flow_conservation][s_node, cut_edge], new_var, -1/v_subgraph.nb_appearance_nodes[v_node] )
        end

        if v_node ∈ v_subgraph.overlapping_nodes
            idx = v_subgraph.idx_overlapping[v_node]
            set_normalized_coefficient(model_master[:overlapping][s_node, v_node, idx], new_var, -1)
            if idx < v_subgraph.nb_appearance_nodes[v_node]
                set_normalized_coefficient(model_master[:overlapping][s_node, v_node, idx+1], new_var, +1)
            else
                set_normalized_coefficient(model_master[:overlapping][s_node, v_node, 1], new_var, +1)
            end
        end

    end
    
    for (i_edge, v_edge) in enumerate(v_subgraph.edges)
        path = submapping.edge_routing[i_edge]
        curr_ve_dem = ve_dem[src(v_edge), dst(v_edge)]
        for i_node in 1:length(path)-1
            src = path[i_node]
            dst = path[i_node+1]
            if dst < src
                src, dst = dst, src
            end
            cost += curr_ve_dem * se_cost[src, dst]
            set_normalized_coefficient(model_master[:edge_capacity][Edge(src, dst)], new_var, curr_ve_dem)
        end
    end
    set_objective_coefficient(model_master, new_var, cost)
    set_normalized_coefficient(model_master[:submapping_selection][v_subgraph], new_var, 1)
end



function add_dumb_columns!(model_master::Model, v_decomposition::OverlappingVirtualDecomposition)
    for v_subgraph in v_decomposition.subgraphs
        dumb_var = @variable(model_master, lower_bound=0., upper_bound=1.)
        set_objective_coefficient(model_master, dumb_var, 10e4)
        set_normalized_coefficient(model_master[:submapping_selection][v_subgraph], dumb_var, 1)
    end
end


function add_single_node_columns!(model_master::Model, columns::Dict{Subgraph, Vector{Column}}, instance::Instance, v_decomposition::OverlappingVirtualDecomposition)
    
    v_g, vn_dem, ve_dem = instance.v_network.graph, instance.v_network.node_demands, instance.v_network.edge_demands
    s_g, s_dir, sn_cap, se_cap, sn_cost, se_cost = instance.s_network.graph, instance.s_network.directed_graph, instance.s_network.node_capacities, instance.s_network.edge_capacities, instance.s_network.node_costs, instance.s_network.edge_costs

    for v_subgraph in v_decomposition.subgraphs
        if isone(length(v_subgraph.nodes))
            v_node = v_subgraph.nodes[1]
            for s_node in vertices(s_g)
                if sn_cap[s_node] >= vn_dem[v_node]
                    placement = zeros(Int, nv(v_g))
                    placement[v_node] = s_node
                    submapping = Mapping(placement, Vector{Vector{Int}}())
                    add_column!(model_master, columns, v_subgraph, submapping, instance)
                end
            end
        end
    end

end




struct DualValues
    submapping_selection::JuMP.Containers.DenseAxisArray
    overlapping::JuMP.Containers.SparseAxisArray
    node_1t1::AbstractArray{Float64}
    edge_capacity::JuMP.Containers.DenseAxisArray
    flow_conservation::JuMP.Containers.DenseAxisArray
    flow_departure::JuMP.Containers.DenseAxisArray
end


function DualValues(model::Model)
    return DualValues(
        dual.(model[:submapping_selection]),
        dual.(model[:overlapping]),
        dual.(model[:node_1t1]),
        dual.(model[:edge_capacity]),
        dual.(model[:flow_conservation]),
        dual.(model[:flow_departure])
    )
end

# doesnt work
function zero_duals(model::Model)
    return DualValues(
        JuMP.Containers.DenseAxisArray(zeros(size(model[:submapping_selection])), axes(model[:submapping_selection])...),
        JuMP.Containers.DenseAxisArray(zeros(size(model[:overlapping])), axes(model[:overlapping])...),
        zeros(size(model[:node_1t1])),
        JuMP.Containers.DenseAxisArray(zeros(size(model[:edge_capacity])), axes(model[:edge_capacity])...),
        JuMP.Containers.DenseAxisArray(zeros(size(model[:flow_conservation])), axes(model[:flow_conservation])...),
        JuMP.Containers.DenseAxisArray(zeros(size(model[:flow_departure])), axes(model[:flow_departure])...)
    )
end


function stabilize_duals!(dual_costs::DualValues, current_dual_costs::DualValues, stabilization_coeff::Float64)

    dual_costs.submapping_selection.data .= dual_costs.submapping_selection.data * stabilization_coeff + current_dual_costs.submapping_selection.data * (1 - stabilization_coeff)
    dual_costs.node_1t1 .= dual_costs.node_1t1 * stabilization_coeff + current_dual_costs.node_1t1 * (1 - stabilization_coeff)
    dual_costs.flow_conservation.data .= dual_costs.flow_conservation.data * stabilization_coeff + current_dual_costs.flow_conservation.data * (1 - stabilization_coeff)
    dual_costs.edge_capacity.data .= dual_costs.edge_capacity.data * stabilization_coeff + current_dual_costs.edge_capacity.data * (1 - stabilization_coeff)
    dual_costs.flow_departure.data .= dual_costs.flow_departure.data * stabilization_coeff + current_dual_costs.flow_departure.data * (1 - stabilization_coeff)


    # 3. Sparse Axis Array (Loop through keys of the OrderedDict to avoid allocations)
    α = stabilization_coeff
    β = 1.0 - α
    
    # We retrieve the underlying dicts
    dict_dest = dual_costs.overlapping.data
    dict_src = current_dual_costs.overlapping.data
    
    for (key, val) in dict_dest
        # Update in-place using the corresponding value from current_dual_costs
        dict_dest[key] = val * α + dict_src[key] * β
    end
end