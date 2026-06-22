
struct VirtualNetwork
    graph::Graph
    name::String
    node_demands::Vector{Int}
    edge_demands::Matrix{Int}
end

struct SubstrateNetwork
    graph::Graph
    name::String
    directed_graph::DiGraph
    node_capacities::Vector{Int}
    node_costs::Vector{Int}
    edge_capacities::Matrix{Int}
    edge_costs::Matrix{Int}
end


struct InstanceVNE
    v_network::VirtualNetwork
    s_network::SubstrateNetwork
end


struct Mapping 
    node_placement::Vector{Int}
    edge_routing::Vector{Vector{Int}}
end