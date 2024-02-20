import Oscar.TropicalSemiringElem

abstract type Strategy end

struct Series <: Strategy end
struct Step <: Strategy end
struct Parallel <: Strategy end

struct MixedCellTracker{strategy<:Strategy}
    pointers::Vector{Int}
    dual_cell_trackers::Vector{DualCell}
end

function mixed_cell_tracker(dual_cell_trackers::Vector{DualCellTracker}, strategy::Strategy)
    check_mixed_cell_tracker_inputs(dual_cell_trackers)
    return MixedCellTracker{strategy}(zeros(length(dual_cell_trackers), 1), dual_cell_trackers)
end

function check_mixed_cell_tracker_inputs(dual_cell_trackers)
    # to do: implement
end

function Base.show(io::IO, T::MixedCellTracker)
    print(io, "Mixed cell tracker with pointers $(T.pointers) and dual cell trackers $(T.dual_cell_trackers)")
end