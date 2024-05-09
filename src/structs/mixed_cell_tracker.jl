import Oscar.TropicalSemiringElem

"""
This struct couples a mixed cell to a mixed path, which allows one to ask for things like breaking points or drifts. This data forms the fundamental building blocks of our tropical homotopy algorithm.

The mixed cell is assumed to arise from a generic lift between the head of the mixed path and the point of interest thereafter. Thus, it always makes sense to ask for the corresponding transverse intersection point, drift, or next breaking point (if it exists).
"""
struct MixedCellTracker
    mixed_path::MixedPath
    mixed_cell::MixedCell
end

function mixed_cell_tracker(mixed_path::MixedPath, mixed_cell::MixedCell)
    return MixedCellTracker(mixed_path, mixed_cell)
end

function Base.show(io::IO, T::MixedCellTracker)
    print(io, "Tracker for a mixed cell along a mixed path with mixed vector $(mixed_vector(T))")
end

function mixed_cell(T::MixedCellTracker)::MixedCell
    return T.mixed_cell
end

function mixed_path(T::MixedCellTracker)::MixedPath
    return T.mixed_path
end

function mixed_vector(T::MixedCellTracker)
    return mixed_vector(mixed_cell(T))
end

function pointers(T::MixedCellTracker)
    return pointers(mixed_path(T))
end