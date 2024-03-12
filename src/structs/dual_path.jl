using Oscar

abstract type DualPathType end

struct DualPathHypersurface <: DualPathType end
struct DualPathLinear <: DualPathType end
struct DualPathInvertedLinear <: DualPathType end


struct DualPath{pathType<:DualPathType, minOrMax<:Union{typeof(min),typeof(max)}}

    nodes::Vector{Vector{Oscar.TropicalSemiringElem{minOrMax}}}

    function DualPath{pathType, minOrMax}(nodes::Vector{Vector{Oscar.TropicalSemiringElem{minOrMax}}}) where {pathType<:DualPathType, minOrMax<:Union{typeof(min),typeof(max)}}
        return new{pathType, minOrMax}(nodes)
    end

end