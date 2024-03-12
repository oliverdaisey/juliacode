
"""
    stable_intersection_point(s::MixedCell, h::MixedPath, pointer::Vector{Int}, fraction::QQFieldElem)

Given a mixed cell `s`, a mixed path `h`, an index representing the pointer to the most recently visited node, and a fraction between the most recently visited node and the next node, return the intersection point of the stable intersection of `s` and `h` at the given pointer and fraction.

Note that the data of a pointer_index with a fraction t uniquely determines a lift.
*"""
function stable_intersection_point(s::MixedCell, h::MixedPath, pointer_index::Int, t::QQFieldElem)
    # write down the matrix of equations determining the intersection point
    
    coefficients = Vector{Vector{Int}}()
    constants = Vector{QQFieldElem}()
    for i in 1:length(s.dual_cells)
        dualCell = s.dual_cells[i]
        # iterate through all 2-element subsets of dualCell.activeSupport
        for subset in subsets(dualCell.activeSupport, 2)
            k, l = subset
            # coefficients of the new equation
            push!(coefficients, dualCell.ambientSupport[k, :] - dualCell.ambientSupport[l, :])
            # constant of the new equation
            # first get the required components of the lifts at this time
            h_k = lift_from_node_and_fraction(h.dualPaths[i], h.pointers[pointer_index][i], t)[k]
            h_l = lift_from_node_and_fraction(h.dualPaths[i], h.pointers[pointer_index][i], t)[l]
            push!(constants, h_l - h_k)
        end
    end

    return solve(matrix(QQ, coefficients), matrix(QQ, hcat(constants)))

end