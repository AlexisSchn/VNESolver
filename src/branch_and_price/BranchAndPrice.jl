module BranchAndPrice

using ..VNESolver
using JuMP, CPLEX, DataStructures, Graphs, Random, Printf


export BranchAndPriceResult, BranchAndPriceParameters
export solve_branch_price

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

include("network_decomposition.jl")
include("branch_and_price.jl")
include("branching/branching.jl")
include("branching/first_try.jl")
include("master_problem.jl")
include("column_generation/root_node.jl")
include("column_generation/tree_node.jl")
include("pricers/milp.jl")
include("pricers/greedy.jl")
include("pricers/greedy_substrate_subgraph.jl")
include("pricers/greedy_branching.jl")




end # module