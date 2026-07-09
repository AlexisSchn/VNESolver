module FlowFormulation

using ..VNESolver
using JuMP, CPLEX, Graphs


struct FlowFormulationResult <: AbstractSolverResult
    vn_name::String
    sn_name::String
    objective_value::Float64
    lower_bound::Float64
    gap::Float64
    nb_nodes::Int
    solving_time::Float64
    status::Symbol
end


function FlowFormulationResult(instance::Instance, model::Model, status::Symbol)
    if status==:Feasible
        return FlowFormulationResult(
            instance.v_network.name, 
            instance.s_network.name, 
            objective_value(model), 
            objective_bound(model),
            relative_gap(model),
            node_count(model), 
            solve_time(model), 
            :Feasible
        )
    elseif status==:Unfeasible
        return FlowFormulationResult(
            instance.v_network.name, 
            instance.s_network.name, 
            Inf, 
            objective_bound(model),
            Inf, 
            node_count,
            solve_time(model), 
            :Unfeasible
        )
    end
end

include("flow_formulation.jl")

export solve_flow_formulation, solve_flow_formulation_linear

end # module