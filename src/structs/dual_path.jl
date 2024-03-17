"""
    DualPath{pathType, minOrMax}

A dual path is a sequence of nodes, each of which is a vector of tropical semiring elements. The nodes are connected by a linear interpolation. The ambient dual support is the support of the dual cell. Construct them with this data, or alternatively use the `dual_path` function.
"""
struct DualPath{pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}

    nodes::Vector{Vector{Oscar.TropicalSemiringElem{minOrMax}}}
    ambientDualSupport::DualSupport{pathType}

    function DualPath{pathType, minOrMax}(nodes::Vector{Vector{Oscar.TropicalSemiringElem{minOrMax}}}, ambientDualSupport::DualSupport{pathType}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}
        return new{pathType, minOrMax}(nodes, ambientDualSupport)
    end

end

"""
    dual_path(nodes::Vector{Vector{Oscar.TropicalSemiringElem{minOrMax}}}, ambientDualSupport::DualSupport{pathType}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}

Create a dual path from the given nodes and ambient dual support.
"""
function dual_path(nodes::Vector{Vector{Oscar.TropicalSemiringElem{minOrMax}}}, ambientDualSupport::DualSupport{pathType}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}
    return DualPath{pathType, minOrMax}(nodes, ambientDualSupport)
end

"""
    dual_path(d::DualCell{pathType, minOrMax}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}

Create a dual path from a dual cell.
"""
function dual_path(d::DualCell{pathType, minOrMax}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}
    return DualPath{pathType, minOrMax}([dual_vector(d)], ambient_dual_support(d))
end

"""
    lift_from_node_and_fraction(h::DualPath, index::Int, t::QQFieldElem)

Get the dual weight from a node of a dual path indexed by `index` and a fraction of the distance to the next node.
"""
function lift_from_node_and_fraction(h::DualPath, index::Int, t::QQFieldElem)

    @assert 1 <= index <= length(h.nodes) "Index out of bounds"
    @assert 0 <= t <= 1 "t out of bounds"

    if index >= length(h.nodes)
        return [x.data for x in h.nodes[index]]
    end

    return [x.data for x in h.nodes[index]] + t * ([x.data for x in h.nodes[index+1]] - [x.data for x in h.nodes[index]])
end

function Base.show(io::IO, h::DualPath)
    print(io, "Dual path with nodes $(h.nodes)")
end

function Base.getindex(h::DualPath, i::Int)
    return h.nodes[i]
end

function nodes(h::DualPath)
    return h.nodes
end