
"""
    next_point_of_interest(T::MixedCellTracker)

Given a mixed cell tracker, returns the next point of interest (either the next node, or the breaking point if it exists), along with the supports that change.

"""
function next_point_of_interest(T::MixedCellTracker)

    # make this code work with the old implementation (fix this after deadline)
    h = T.mixed_path
    pointer_index = 1
    fraction = QQFieldElem(0)

    lengths = length.(ambient_support.(dual_cells(mixed_cell(T))))

    # deal with the silly case that the mixed path only has one node left
    if pointer_index == length(pointers(h))
        return nothing, [] # this means that this mixed cell tracker is done
    end


    dual_path_pointers = pointers(h)[pointer_index]
    n = length(dual_path_pointers)

    # get the lift at this time
    lift = vcat([lift_from_node_and_fraction(dual_paths(h)[i], dual_path_pointers[i], fraction) for i in 1:n]...)

    # get direction path is travelling in dual space
    next_dual_path_pointers = pointers(h)[pointer_index+1]
    direction = vcat([nodes(dual_paths(h)[i])[next_dual_path_pointers[i]] for i in 1:n]...) ./ vcat([nodes(dual_paths(h)[i])[dual_path_pointers[i]] for i in 1:n]...)

    C_s = mixed_cell_cone(s)

    # starting at `lift` and moving in the direction `direction`, when do we hit facet of C_s?
    # we want to find the smallest t such that lift + t*direction is on the facet of C_s
    t = ray_intersects_cone(mixed_cell_cone_to_polyhedron(C_s), lift, QQ.(direction))


    if t > 1 || isnothing(t)
        # this means that the breaking point does not exist or takes us away from the path, so we should return the next node (no supports change)
        return partition_vector(lengths, h[2]), []
    end

    # otherwise we are in the case that the breaking point exists, work out which supports change
    facet_point = lift + t*QQ.(direction)

    # get supports that change
    support_indices = []
    for i in 1:length(coefficients(C_s))
        if iszero(facet_point .* coefficients(C_s)[i])
            # intersects this facet, supports that change are indices of nonzero entries of coefficients of C_s
            push!(support_indices, findall(!iszero, coefficients(C_s)[i]))
        end
    end

    # make sure support_indices is unique
    support_indices = unique(vcat(support_indices...))

    return partition_vector(lengths, TT.(facet_point)), support_indices

end

function partition_vector(lengths::Vector{Int}, v::Vector{Oscar.TropicalSemiringElem{minOrMax}})::Vector{Vector{Oscar.TropicalSemiringElem{minOrMax}}} where {minOrMax<:Union{typeof(min),typeof(max)}}
    # first vector should have length lengths[1], second vector should have length lengths[2], etc.

    partitionedVector = Vector{Oscar.TropicalSemiringElem{minOrMax}}[]
    startIndex = 1
    for length in lengths
        endIndex = startIndex + length - 1
        push!(partitionedVector, v[startIndex:endIndex])
        startIndex = endIndex + 1
    end
    return partitionedVector

end