
function tropical_drift(T::MixedCellTracker)

    # to do: infer a correct epsilon analytically
    epsilon = QQFieldElem(1//100)
    initialPoint = stable_intersection_point(T.mixed_cell)
    dualDirection = mixed_path(T)[2] ./ mixed_path(T)[1]

    lengths = [length(dual_weight(x)) for x in dual_cells(mixed_cell(T))]

    # now I want to split dualDirection into its components given by lengths
    # so if lengths = [2, 6] and dualDirection = [1, 2, 3, 4, 5, 6, 7, 8] then I want to get [[1, 2], [3, 4, 5, 6, 7, 8]]
    dualDirectionVectors = [dualDirection[i:i+lengths[j]-1] for (i, j) in zip([1; cumsum(lengths)[1:end-1] .+ 1], 1:length(lengths))]

    # build a new mixed cell using a lift a small distance along this dual direction (small enough to remain in mixed cell cone)
    newMixedCell = deepcopy(mixed_cell(T))
    for i in 1:length(dual_cells(newMixedCell))
        newMixedCell.dual_cells[i].dualWeight = newMixedCell.dual_cells[i].dualWeight .* TT.(epsilon * [x.data for x in dualDirectionVectors[i]])
    end

    finalPoint = stable_intersection_point(newMixedCell)

    return finalPoint - initialPoint
end