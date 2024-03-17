using Oscar

"""
    dual_cells(S::DualSupport{<:DualType}, c::Vector{Oscar.TropicalSemiringElem})

Create a vector of all the dual cells from a support (of exponent vectors) and a lift.

"""
function dual_cells(S::DualSupport{<:DualType}, c::Vector{<:Oscar.TropicalSemiringElem})

    dualSubdivision = subdivision_of_points_workaround(S.points, c)
    polyhedralComplex = polyhedral_complex(Oscar.pm_object(dualSubdivision).POLYHEDRAL_COMPLEX)

    dualCells = []
    for i in tropical_codim(S):dim(polyhedralComplex)
        for p in polyhedra_of_dim(polyhedralComplex, i)
            s = [Int.(v) for v in vertices(p)]
            if is_dual_cell_candidate(S, s)
                push!(dualCells, dual_cell(S, s, c))
            end
        end
    end

    return dualCells
end