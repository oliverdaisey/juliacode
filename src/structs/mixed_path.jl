include("../structs/dual_path.jl")

struct MixedPath

    pointers::Vector{Vector{Int}}
    dualPaths::Vector{<:DualPath}

    function MixedPath(pointers::Vector{Vector{Int}}, dualPaths::Vector{<:DualPath})
        return new(pointers, dualPaths)
    end

end

function mixed_path_in_series(dualPaths::Vector{<:DualPath})
    
    # create the list of pointers
    lengths = [length(dualPath.nodes) for dualPath in dualPaths]

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