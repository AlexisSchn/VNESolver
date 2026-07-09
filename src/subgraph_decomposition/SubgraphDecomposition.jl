module SubgraphDecomposition

using ..VNESolver
using JuMP, CPLEX, DataStructures, Graphs, Random, Printf


export solve_subgraph_decomposition, SubgraphDecompositionResult, SubgraphDecompositionParameters


struct SubgraphDecompositionResult <: AbstractSolverResult
    vn_name::String
    sn_name::String
    rmp_value::Float64
    lg_bound::Float64
    gap::Float64
    nb_columns::Int
    nb_iter::Int
    solving_time::Float64
end

struct SubgraphDecompositionParameters <: AbstractSolverParameters
    time_max::Float64
    nb_iter_max::Int
    nb_columns_max::Int
    gap_min::Float64
end

function SubgraphDecompositionParameters()
    return SubgraphDecompositionParameters(
        500.,
        500,
        5000,
        0.01
    )
end

include("master_problem.jl")
include("pricer.jl")
include("pricer_greedy.jl")
include("column_generation.jl")




end # module