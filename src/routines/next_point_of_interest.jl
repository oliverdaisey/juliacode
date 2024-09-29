
"""
    next_point_of_interest(T::MixedCellTracker)

Given a mixed cell tracker, returns the next point of interest (either the next node, or the breaking point if it exists, or `nothing` otherwise), along with the supports that change.

Has a side effect: Changes the dual weights in the dual cells of the mixed cell `T` is tracking, to the point of interest (if it exists). Also updates the path of the mixed cell tracker when required. (Basically, this function is doing too much, for now).

Do not call this function unless you are actually intending to traverse the mixed path!
"""
function next_point_of_interest(T::MixedCellTracker)

    # make this code work with the old implementation (fix this after deadline)
    h = T.mixed_path
    pointer_index = 1
    fraction = QQFieldElem(0)
    lengths = length.(ambient_support.(dual_cells(mixed_cell(T))))

    # deal with case that the mixed path only has one node left
    if length(pointers(h)) == 1
        return nothing, [] # this means that this mixed cell tracker is done
    end


    dual_path_pointers = pointers(h)[1]
    n = length(dual_path_pointers)

    # get the lift at this time
    lift = QQ.(vcat(to_vector.(mixed_vector(T))...))

    # get direction path is travelling in dual space
    next_dual_path_pointers = pointers(h)[2]

    C_s = mixed_cell_cone(s)

    # starting at `lift` and moving in the direction `direction`, when do we hit facet of C_s?
    # we want to find the smallest t such that lift + t*direction is on the facet of C_s
    t = ray_intersects_cone(mixed_cell_cone_to_polyhedron(C_s), lift, QQ.(direction(h, 1)))

    if isnothing(t) || t > 1 
        # this means that the breaking point does not exist or takes us away from the path, so we should return the next node (no supports change)
        partitionedVector = to_vector.(h[2])
        update_dual_weights!(T, partitionedVector)
        update_path(T)
        return partitionedVector, []
    end

    # otherwise we are in the case that the breaking point exists, work out which supports change
    facet_point = lift + t*QQ.(direction(h,1))

    # get supports that change
    support_indices = []
    for i in 1:length(inequalities(C_s))
        if iszero(facet_point .* inequalities(C_s)[i])
            # intersects this facet, supports that change are indices of nonzero entries of coefficients of C_s
            push!(support_indices, findall(iszero, inequalities(C_s)[i]))
        end
    end

    # make sure support_indices is unique
    support_indices = unique(vcat(support_indices...))

    # update dual weights
    update_dual_weights!(T, partition_vector(lengths, TT.(facet_point)))
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

function update_dual_weights!(T::MixedCellTracker, newDualWeights::Vector{Vector{Oscar.TropicalSemiringElem{minOrMax}}}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    s = mixed_cell(T)
    for i in 1:length(dual_cells(s))
        update!(s.dual_cells[i].dualWeight, newDualWeights[i])
    end
end

function update_path(T::MixedCellTracker)
    # pop the first pointer from the path
    mixed_path(T).pointers = mixed_path(T).pointers[2:end]
end