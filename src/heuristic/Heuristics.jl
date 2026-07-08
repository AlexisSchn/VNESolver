module Heuristics

using ..VNESolver
using DataStructures, Graphs

include("base.jl")
include("greedy.jl")
include("local_search.jl")

export solve_greedy, local_search

end # module