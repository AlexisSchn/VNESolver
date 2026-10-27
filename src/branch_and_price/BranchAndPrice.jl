module BranchAndPrice

using ..VNESolver
using JuMP, CPLEX, DataStructures, Graphs, Random, Printf



export BranchAndPriceResult, BranchAndPriceParameters
export solve_branch_price

include("types.jl")  


include("network_decomposition.jl")
include("branching/branching.jl")
include("branch_and_price.jl")
include("master_problem.jl")
include("column_generation/root_node.jl")
include("column_generation/tree_node.jl")
include("column_generation/potential_solution.jl")
include("pricers/utils_pricers.jl")
include("pricers/exact_routing.jl")
include("pricers/milp.jl")
include("pricers/greedy.jl")
include("pricers/exact.jl")

end # module