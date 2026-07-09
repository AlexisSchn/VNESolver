module FlowFormulation

using ..VNESolver
using JuMP, CPLEX, Graphs

include("flow_formulation.jl")

struct FlowFormulationResult <: AbstractSolverResult
    vn_name::String
    sn_name::String
    objective_value::Float64
    lower_bound::Float64
    gap::Float64
    solving_time::Float64
    status::Symbol
end

function FlowFormulationResult(instance::Instance)
    return FlowFormulationResult("", "", 30., 30., 3., 35., :Feasible)
end

export solve_flow_formulation, solve_flow_formulation_linear

end # module