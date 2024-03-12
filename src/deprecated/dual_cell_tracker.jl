import Oscar.TropicalSemiringElem

struct DualCellTracker
    dual_cell::DualCell
    nodes::Vector{Vector{TropicalSemiringElem}}

    function DualCellTracker(dual_cell::DualCell, nodes::Vector{Vector{TropicalSemiringElem}})
        return new(dual_cell, nodes)
    end
end