
"""
    Support{dualType}(points::Matrix{Int})

A support for a dual object, with rows as points.

"""
mutable struct Support{dualType<:DualType}
    points::Matrix{Int}

    function Support{dualType}(points::Matrix{Int}) where {dualType<:DualType}
        return new{dualType}(points)
    end
end

function points(s::Support)
    return s.points
end

function tropical_ambient_dim(s::Support)
    return size(s.points, 2)
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

function Base.length(s::Support)
    return size(s.points, 1)
end

function Base.iterate(s::Support)
    return iterate(s[i] for i in 1:length(s))
end