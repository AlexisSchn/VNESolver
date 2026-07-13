module VNESolver

using Graphs, JSON
#using GraphRecipes, Plots

export Instance, VirtualNetwork, SubstrateNetwork, Mapping, AbstractSolverResult, AbstractSolverParameters
export get_instance_from_folder, read_substrate, read_virtual, read_virtuals_folder, read_substrates_folder
export get_mapping_cost
export visu_graph, visu_partitioning

# Local
include("core/types.jl")
include("core/tools.jl")
include("core/io.jl")
#include("core/visualization.jl")

abstract type AbstractSolverResult end
abstract type AbstractSolverParameters end

# Submodules
include("compact/FlowFormulation.jl")
using .FlowFormulation
export solve_flow_formulation, solve_flow_formulation_linear

include("heuristic/Heuristics.jl")
using .Heuristics
export solve_greedy, local_search

include("path_formulation/PathFormulation.jl")
using .PathFormulation
export solve_path_formulation

include("subgraph_decomposition/SubgraphDecomposition.jl")
using .SubgraphDecomposition
export solve_subgraph_decomposition, solve_subgraph_decomposition_better

include("graph_partitioning/GraphPartitioning.jl")
using .GraphPartitioning
export partition_kahip, partition_metis

end # module