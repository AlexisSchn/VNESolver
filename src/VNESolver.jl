module VNESolver

using Random

using JuMP
using Graphs
using CPLEX
using DataStructures
import JSON

export Instance, VirtualNetwork, SubstrateNetwork, Mapping 
export get_instance_from_folder, read_substrate, read_virtual 
export get_mapping_cost

# Local
include("core/types.jl")
include("core/tools.jl")
include("core/io.jl")


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
export solve_subgraph_decomposition

end # module