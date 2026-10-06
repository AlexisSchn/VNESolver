


struct Branching
    placement::Dict{Int, Int}
    zoning::Dict{Int, Vector{Int}}
end


struct TreeNode
    id::Int
    depth::Int
    branching::Branching
    lower_bound::Float64
end