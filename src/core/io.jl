# src/core/io.jl
# Utilities for parsing JSON instance files into core graph structures.

function read_substrate(json_graph)
    
    g = Graph()
    dir_g = DiGraph()
    name = json_graph["name"]

    num_nodes = length(json_graph["nodes"])
    node_capacities = Vector{Int}(undef, num_nodes)
    node_costs = Vector{Int}(undef, num_nodes)
    for (i, node) in enumerate(json_graph["nodes"])
        add_vertex!(g)
        add_vertex!(dir_g)
        node_capacities[i] = node["cap"]
        node_costs[i] = node["cost"]
    end

    edge_capacities = zeros(Int, num_nodes, num_nodes)
    edge_costs = zeros(Int, num_nodes, num_nodes)
    edge_ids = zeros(Int, num_nodes, num_nodes)
    idx = 1
    for edge in json_graph["edges"]
        src = edge["source"]
        dst = edge["target"]
        if src == dst
            error("JSON Parser Error: Loop in the substrate network on node $src.")
        end
        if src > dst 
            @warn "Some undirected nodes are not well ordered on edge ($src, $dst). Be careful please!"    
            temp = src
            src = dst 
            dst = temp
        end

        add_edge!(g, src, dst)
        add_edge!(dir_g, src, dst)
        add_edge!(dir_g, dst, src)
        edge_capacities[src, dst] = edge["cap"]
        edge_capacities[dst, src] = edge["cap"]
        edge_costs[src, dst] = edge["cost"]
        edge_costs[dst, src] = edge["cost"]
        edge_ids[src, dst] = idx
        edge_ids[dst, src] = idx
        idx += 1
    end

    
    if !is_connected(g)
        error("Input Error: The substrate network graph is disconnected! Parts of the graph:$(connected_components(g))")
    end

    return SubstrateNetwork(g, name, dir_g, node_capacities, node_costs, edge_capacities, edge_costs, edge_ids)
end


function read_virtual(json_graph)

    g = Graph()
    name = json_graph["name"]

    num_nodes = length(json_graph["nodes"])
    node_demands = Vector{Int}(undef, num_nodes)
    for (i, node) in enumerate(json_graph["nodes"])
        add_vertex!(g)
        node_demands[i] = node["dem"]
    end

    edge_demands = zeros(Int, num_nodes, num_nodes)
    edge_ids = zeros(Int, num_nodes, num_nodes)
    idx = 1
    for edge in json_graph["edges"]
        src = edge["source"]
        dst = edge["target"]
        if src == dst
            error("JSON Parser Error: Loop in the virtual network on node $src.")
        end
        if src > dst 
            @warn "Some undirected nodes are not well ordered on edge ($src, $dst). Be careful please!"    
            temp = src
            src = dst 
            dst = temp
        end

        add_edge!(g, src, dst)
        edge_demands[src, dst] = edge["dem"]
        edge_ids[src, dst] = idx
        idx += 1
    end
    
    
    if !is_connected(g)
        error("Input Error: The virtual network graph is disconnected! Parts of the graph:$(connected_components(g))")
    end

    return VirtualNetwork(g, name, node_demands, edge_demands, edge_ids)
end



function get_instance_from_folder(folder_path::String)

    virtual_network::Union{Nothing, VirtualNetwork} = nothing
    substrate_network::Union{Nothing, SubstrateNetwork} = nothing    
    
    for filename in readdir(folder_path; join=true)
        isdir(filename) && continue
        !endswith(filename, ".json") && continue

        json_graph = JSON.parsefile(filename)

        if json_graph["type"] == "virtual"
            virtual_network = read_virtual(json_graph)
        elseif json_graph["type"] == "substrate"
            substrate_network = read_substrate(json_graph)
        end
    end


    # --- graph validation ---
    if isnothing(virtual_network) || isnothing(substrate_network)
        error("Could not construct Instance: Missing a virtual or substrate network file in '$folder_path'")
    end
    
    return Instance(virtual_network, substrate_network)
end


function read_virtuals_folder(folderpath)
    vns = VirtualNetwork[]
    for filename in readdir(folderpath; join=true)
        json_graph = JSON.parsefile(filename)
        vn = read_virtual(json_graph)
        push!(vns, vn)
    end
    return vns
end


function read_substrates_folder(folderpath)
    sns = SubstrateNetwork[]
    for filename in readdir(folderpath; join=true)
        json_graph = JSON.parsefile(filename)
        sn = read_substrate(json_graph)
        push!(sns, sn)
    end
    return sns
end

