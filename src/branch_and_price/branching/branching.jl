


struct Branching
    placement_mandatory::Dict{Int, Int}
    placement_allowed::Dict{Int, Vector{Int}}
    routing_mandatory::Dict{Edge, Edge}
    routing_forbidden::Dict{Edge, Vector{Edge}}
end

Branching() = Branching(
    Dict{Int, Int}(),
    Dict{Int, Vector{Int}}(),
    Dict{Edge, Edge}(),
    Dict{Edge, Vector{Edge}}()
)


struct TreeNode
    id::Int
    depth::Int
    branching::Branching
    lower_bound::Float64
end