include("dual_type.jl")

@doc raw"""
    Support{dualType}(points::Matrix{Int})

    A support for a dual object, with rows as points.

"""
mutable struct Support{dualType<:DualType}
    points::Matrix{Int}

    function Support{dualType}(points::Matrix{Int}) where {dualType<:DualType}
        return new{dualType}(points)
    end
end

function Base.getindex(s::Support, i::Int)
    return s.points[i, :]
end

function Base.getindex(s::Support, i::Int, j::Int)
    return s.points[i, j]
end

function Base.getindex(s::Support, I::Vector{Int})
    return s.points[I, :]
end

function Base.getindex(s::Support, i::Int, ::Colon)
    return s.points[i, :]
end

function Base.size(s::Support, i::Int)
    return size(s.points, i)
end