struct MixedCell{minOrMax<:Union{typeof(min),typeof(max)}}
    dual_cells::Vector{DualCell}
end

"""
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

    @assert length(dualCells) > 0 "Mixed cell must contain at least one dual cell"
    @assert length(unique(convention.(dualCells))) == 1 "All dual cells must use the same min/max convention"
    @assert length(unique(ambient_dim.(dualCells))) == 1 "All dual cells must have the same ambient dimension"

    n = ambient_dim(dualCells[1])

    @assert sum(tropical_codim.(dualCells)) + tropical_lineality_dim(dualCells) == n "Dual cells must have complementary dimensions"
end

function ambient_support(s::MixedCell)
    return vcat([dualCell.ambientDualSupport for dualCell in s.dual_cells]...)
end

function Base.show(io::IO, s::MixedCell)
    print(io, "Mixed cell with dual cells $(s.dual_cells)")
end

function Base.copy(s::MixedCell)
    return MixedCell{typeof(min)}(copy(s.dual_cells))
end

function dual_cells(s::MixedCell)
    return s.dual_cells
end

function mixed_vector(s::MixedCell)
    return dual_vector.(dual_cells(s))
end

function tropical_lineality_dim(dualCells::Vector{<: DualCell})
    hulls = convex_hull.(points.(ambient_support.(dualCells)))
    return codim(sum(hulls))
end

function tropical_lineality_space(s::MixedCell)
    sumHulls = sum(convex_hull.(points.(ambient_support.(dual_cells(s)))))
    
    return affine_equation_matrix(affine_hull(sumHulls))[:,2:end]

end