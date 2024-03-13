import Oscar.TropicalSemiringElem
include("mixed_cell.jl")
include("mixed_path.jl")

struct MixedCellTracker
    mixed_path::MixedPath
    mixed_cell::MixedCell
end

function mixed_cell_tracker(mixed_path::MixedPath, mixed_cell::MixedCell)
    @assert ambient_support(mixed_path) == ambient_support(mixed_cell) "Ambient supports do not match"
    return MixedCellTracker(mixed_path, mixed_cell)
end

function check_mixed_cell_tracker_inputs(dual_cell_trackers)
    # to do: implement
end

function Base.show(io::IO, T::MixedCellTracker)
    print(io, "Mixed cell tracker with pointers $(T.pointers) and dual cell trackers $(T.dual_cell_trackers)")
end