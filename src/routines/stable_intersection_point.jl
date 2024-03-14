
"""
    stable_intersection_point(s::MixedCell)

Returns the tropical intersection point dual to the mixed cell `s`.
"""
function stable_intersection_point(s::MixedCell)

    # we will construct the matrix of equations that determines the intersection point
    exponents = Vector{Vector{Int}}()
    coefficients = Vector{Oscar.TropicalSemiringElem}()
    for i in 1:length(s.dual_cells)
        dualCell = s.dual_cells[i]
        # iterate through all 2-element subsets of dualCell.activeSupport
        for subset in subsets(dualCell.activeSupport, 2)
            k, l = subset
            # coefficients of the new equation
            push!(exponents, dualCell.ambientSupport[k, :] - dualCell.ambientSupport[l, :])
            # constant of the new equation
            push!(coefficients, dualCell.dualVector[l] / dualCell.dualVector[k])
        end
    end

    M = matrix(QQ, exponents)
    @assert rank(M) == ncols(M) "The intersection is not transverse"

    solution = solve(M, matrix(QQ, hcat(coefficients)), side=:right)

    return solution

end