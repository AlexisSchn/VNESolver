module OverlappingDecomposition

using ..VNESolver
using JuMP, CPLEX, DataStructures, Graphs, Random, Printf


export solve_overlapping_decomposition, solve_overlapping_decomposition_simple
export OverlappingDecompositionResult, OverlappingDecompositionParameters

struct OverlappingDecompositionResult <: AbstractSolverResult
    vn_name::String
    sn_name::String
    rmp_value::Float64
    lg_bound::Float64
    gap::Float64
    nb_columns::Int
    nb_iter::Int
    solving_time::Float64
end

struct OverlappingDecompositionParameters <: AbstractSolverParameters
    time_max::Float64
    nb_iter_max::Int
    nb_columns_max::Int
    gap_min::Float64
    stab_coeff::Float64
end

function OverlappingDecompositionParameters()
    return OverlappingDecompositionParameters(
        500.,
        500,
        5000,
        0.05,
        0.
    )
end

include("network_decomposition.jl")
include("master_problem.jl")
include("column_generation.jl")
include("pricers/milp.jl")
include("pricers/greedy.jl")
include("pricers/greedy_substrate_subgraph.jl")




end # module