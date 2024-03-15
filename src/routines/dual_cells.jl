using Oscar

"""
    dual_cells(S::Support{<:DualType}, c::Vector{Oscar.TropicalSemiringElem})

Create a vector of all the maximal dual cells from a support (of exponent vectors) and a lift.

TODO: Return *all* dual cells, not just the maximal ones.
"""
function dual_cells(S::Support{<:DualType}, c::Vector{Oscar.TropicalSemiringElem})

    dual_subdivision = subdivision_of_points(S.points, c)

    return [DualCell(S, active_support, c) for active_support in maximal_cells(dual_subdivision)]
end