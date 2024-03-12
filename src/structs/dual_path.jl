using Oscar

abstract type DualPathType end

struct DualPathHypersurface <: DualPathType end
struct DualPathLinear <: DualPathType end
struct DualPathInvertedLinear <: DualPathType end


struct DualPath{pathType<:DualPathType, minOrMax<:Union{typeof(min),typeof(max)}}

    nodes::Vector{Vector{QQFieldElem}}

    function DualPath{pathType, minOrMax}(nodes::Vector{Vector{QQFieldElem}}) where {pathType<:DualPathType, minOrMax<:Union{typeof(min),typeof(max)}}
        return new{pathType, minOrMax}(nodes)
    end

end

function lift_from_node_and_fraction(h::DualPath, index::Int, t::QQFieldElem)

    @assert 1 <= index <= length(h.nodes) "Index out of bounds"
    @assert 0 <= t <= 1 "t out of bounds"

    if index >= length(h.nodes)
        return h.nodes[index]
    end

    return h.nodes[index] + t * (h.nodes[index+1] - h.nodes[index])
end