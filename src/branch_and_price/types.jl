

# general solver related stuff

struct BranchAndPriceResult <: AbstractSolverResult
    vn_name::String
    sn_name::String
    rmp_value::Float64
    lg_bound::Float64
    gap::Float64
    nb_columns::Int
    nb_iter::Int
    solving_time::Float64
end

struct BranchAndPriceParameters <: AbstractSolverParameters
    time_max::Float64
    nb_iter_max::Int
    nb_columns_max::Int
    gap_min::Float64
    stab_coeff::Float64
end



function BranchAndPriceParameters()
    return BranchAndPriceParameters(
        500.,
        500,
        5000,
        0.05,
        0.
    )
end



# === decomposition

struct Subgraph
    nodes::Vector{Int}
    edges::Vector{Edge}
    cut_edges_src::Vector{Vector{Edge}}
    cut_edges_dst::Vector{Vector{Edge}}
    idx_of_nodes::Vector{Int}
    overlapping_nodes::Vector{Int}
    idx_overlapping::Vector{Int}
    nb_appearance_nodes::Vector{Int}
end

struct OverlappingVirtualDecomposition
    subgraphs::Vector{Subgraph}
    cut_edges::Vector{Edge}
    overlapping_nodes::Vector{Int}
    nb_appearance_nodes::Vector{Int}
end


# solver related

struct Column
    variable::VariableRef
    submapping::Mapping
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



# branch and price related stuff

struct Branching
    placement_mandatory::Dict{Int, Int}
    placement_allowed::Dict{Int, Vector{Int}}
    routing_cut_edges_mandatory::Dict{Edge, Vector{Edge}}
    routing_cut_edges_forbidden::Dict{Edge, Vector{Edge}}
    routing_internal_edges_mandatory::Dict{Edge, Vector{Edge}}
    routing_internal_edges_forbidden::Dict{Edge, Vector{Edge}}
end

Branching() = Branching(
    Dict{Int, Int}(),
    Dict{Int, Vector{Int}}(),
    Dict{Edge, Vector{Edge}}(),
    Dict{Edge, Vector{Edge}}(),
    Dict{Edge, Vector{Edge}}(),
    Dict{Edge, Vector{Edge}}()
)


struct TreeNode
    id::Int
    depth::Int
    branching::Branching
    lower_bound::Float64
end


mutable struct BranchAndPriceState
    instance::Instance
    v_decomposition::OverlappingVirtualDecomposition
    model_master::Model    
    columns
    pricers
    best_solution::Real
    best_bound::Float64
    nb_nodes::Int
    nb_nodes_remaining::Int
    tree_queue::PriorityQueue{TreeNode, Float64}
end



function compare_integer_solution(bnp_state::BranchAndPriceState, new_solution::Int)
    if new_solution < bnp_state.best_solution
        bnp_state.best_solution = new_solution
    end
end




function BranchAndPriceState(instance::Instance, v_decomposition::OverlappingVirtualDecomposition, model::Model, columns, pricers)
    return BranchAndPriceState(
        instance,
        v_decomposition,
        model,
        columns,
        pricers,
        Inf, 
        0.,
        0,
        0,
        PriorityQueue{TreeNode, Float64}()
    )
end
