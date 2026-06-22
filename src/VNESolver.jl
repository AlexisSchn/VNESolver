module VNESolver

using Random

using JuMP
using Graphs
using CPLEX

import JSON


include("core/types.jl")
include("core/tools.jl")
include("core/io.jl")

include("compact/flow_formulation.jl")

include("heuristic/tools.jl")
include("heuristic/greedy.jl")
include("heuristic/local_search.jl")

export Instance, VirtualNetwork, SubstrateNetwork, Mapping 
export get_instance_from_folder, solve_flow_formulation, solve_greedy, local_search
export get_mapping_cost


end # module