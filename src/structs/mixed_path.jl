
abstract type Strategy end

struct Series <: Strategy end
struct Step <: Strategy end
struct Parallel <: Strategy end

mutable struct MixedPath

    pointers::Vector{Vector{Int}}
    dualPaths::Vector{<:DualPath}

    function MixedPath(pointers::Vector{Vector{Int}}, dualPaths::Vector{<:DualPath})
        return new(pointers, dualPaths)
    end

end

function mixed_path_in_series(dualPaths::Vector{<:DualPath})
    
    # create the list of pointers
    lengths = [length(nodes(dualPath)) for dualPath in dualPaths]

    pointers = Vector{Vector{Int}}([[1 for i in 1:length(dualPaths)]])
    
    # initial pointer is [1, 1, 1, ...]
    # final pointer is [length(dualPaths[1].nodes), length(dualPaths[2].nodes), ...]
    # pointers should increase in the first entry until it reaches the length of the first dual path
    # then the second entry should increase until it reaches the length of the second dual path
    # and so on

    for i in 1:length(dualPaths)
        while pointers[end][i] < lengths[i]
            new_pointer = copy(pointers[end])
            new_pointer[i] += 1
            push!(pointers, new_pointer)
        end
    end
    

    return MixedPath(pointers, dualPaths)
end


function lift_from_node_and_fraction(h::MixedPath, index::Int, t::QQFieldElem)

    @assert 1 <= index <= length(h.nodes) "Index out of bounds"
    @assert 0 <= t <= 1 "t out of bounds"
    @assert index != length(h.nodes) && t>0 "Cannot work out lift from the last node"

    return h.nodes[index] + t * (h.nodes[index+1] - h.nodes[index])
end

function ambient_support(h::MixedPath)
    return vcat([dualPath.ambientDualSupport for dualPath in h.dualPaths]...)
end

function Base.show(io::IO, h::MixedPath)
    print(io, "Mixed path with paths of types $(join([type(dualPath) for dualPath in h.dualPaths], ", "))")
end

function Base.getindex(h::MixedPath, i::Int)
    return vcat([h.dualPaths[j][h.pointers[i][j]] for j in 1:length(h.dualPaths)]...)
end

function pointers(h::MixedPath)
    return h.pointers
end

function dual_paths(h::MixedPath)
    return h.dualPaths
end