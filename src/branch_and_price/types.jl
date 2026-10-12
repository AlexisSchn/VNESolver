

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
