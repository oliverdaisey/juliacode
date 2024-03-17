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
    check_mixed_cell_tracker_inputs(mixed_path, mixed_cell)
    return MixedCellTracker(mixed_path, mixed_cell)
end

function check_mixed_cell_tracker_inputs(mixed_path::MixedPath, mixed_cell::MixedCell)
    @assert ambient_support(mixed_path) == ambient_support(mixed_cell) "Ambient supports do not match"
end

function Base.show(io::IO, T::MixedCellTracker)
    print(io, "Tracker for a mixed cell along a mixed path")
end