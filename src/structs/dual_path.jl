"""
    DualPath{pathType, minOrMax}

A dual path is a sequence of nodes, each of which is a vector of tropical semiring elements. The nodes are connected by a linear interpolation. Construct them with this data, or alternatively use the `dual_path` function.
"""
struct DualPath{pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}

    nodes::Vector{DualWeight{pathType, minOrMax}}

    function DualPath{pathType, minOrMax}(nodes::Vector{DualWeight{pathType, minOrMax}}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}
        return new{pathType, minOrMax}(nodes)
    end

end

function dual_path(nodes::Vector{DualWeight{pathType, minOrMax}}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}
    return DualPath{pathType, minOrMax}(nodes)
end

"""
    dual_path(d::DualCell{pathType, minOrMax}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}

Create a singleton dual path from a single dual weight. Use this when you don't want to move one of your dual cells in the homotopy.
"""
function dual_path(w::DualWeight{pathType, minOrMax}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}
    return DualPath{pathType, minOrMax}([w])
end

"""
    dual_path(d::DualCell{pathType, minOrMax}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}

Create a singleton dual path from a dual cell. Use this when you don't want to move one of your dual cells in the homotopy.
"""
function dual_path(d::DualCell{pathType, minOrMax}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}
    return DualPath{pathType, minOrMax}([dual_weight(d)])
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

function Base.getindex(h::DualPath, i::Int)
    return h.nodes[i]
end

function nodes(h::DualPath)
    return h.nodes
end

function type(h::DualPath{pathType, minOrMax}) where {pathType<:DualType, minOrMax<:Union{typeof(min),typeof(max)}}
    return pathType
end

"""
    direction(h::DualPath, index::Int)

Get the direction of the dual path after the given index of a node.
"""
function direction(h::DualPath, index::Int)
    
        @assert 1 <= index <= length(h.nodes) "Index out of bounds"
    
        if index >= length(h.nodes)
            return [0 for x in nodes(h)[1]]
        end
    
        return direction(h[index], h[index+1])
end

function Base.show(io::IO, h::DualPath)
    print(io, "$(type(h)) dual path with $(length(nodes(h))) node(s)")
end