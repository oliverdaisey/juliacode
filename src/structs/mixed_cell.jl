struct MixedCell{minOrMax<:Union{typeof(min),typeof(max)}}
    dual_cells::Vector{DualCell}
end

function mixed_cell(dualCells::Vector{<: DualCell})
    return MixedCell{typeof(min)}(dualCells)
end
