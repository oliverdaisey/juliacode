struct MixedCell{minOrMax<:Union{typeof(min),typeof(max)}}
    dual_cells::Vector{DualCell}
end

r"""
    mixed_cell(dualCells::Vector{<: DualCell})

    Create a mixed cell from the given dual cells.

    # Arguments
    - `dualCells::Vector{<: DualCell}`: The dual cells to be combined into a mixed cell.

    # Returns
    A mixed cell from the given dual cells.

"""
function mixed_cell(dualCells::Vector{<: DualCell})
    check_mixed_cell_inputs(dualCells)
    return MixedCell{typeof(min)}(dualCells)
end

function check_mixed_cell_inputs(dualCells::Vector{<: DualCell})
    # check that the dual cells are using the same convention
    cellTypes = [typeof(dualCell) for dualCell in dualCells]
    @assert all(cellTypes .== cellTypes[2]) "All dual cells must use the same min/max convention"

    @assert length(dualCells) > 0 "Mixed cell must contain at least one dual cell"

    @assert all([size(dualCell.ambientSupport, 1) == size(dualCells[1].ambientSupport, 1) for dualCell in dualCells]) "All dual cells must have the same ambient dimension"

    d = size(dualCells[1].ambientSupport, 1)

    # check that the dual cells are of complementary dimension
    cellDims = [codim(dualCell) for dualCell in dualCells]
    @assert sum(cellDims) == d "Dual cells must have complementary dimensions"
end