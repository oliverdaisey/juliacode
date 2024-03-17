
function tropical_drift(T::MixedCellTracker)

    # to do: infer a correct epsilon analytically
    epsilon = QQFieldElem(1//100)
    initial_point = stable_intersection_point(T.mixed_cell)
    dual_direction = T.mixed_path[2] ./ T.mixed_path[1]

    lengths = [length(x.dualVector) for x in T.mixed_cell.dual_cells]

    # now I want to split dual_direction into its components given by lengths
    # so if lengths = [2, 6] and dual_direction = [1, 2, 3, 4, 5, 6, 7, 8] then I want to get [[1, 2], [3, 4, 5, 6, 7, 8]]
    dual_direction_vectors = [dual_direction[i:i+lengths[j]-1] for (i, j) in zip([1; cumsum(lengths)[1:end-1] .+ 1], 1:length(lengths))]

    # build a new mixed cell using a lift a small distance along this dual direction (small enough to remain in mixed cell cone)
    new_mixed_cell = mixed_cell(T.mixed_cell.dual_cells)
    for i in 1:length(new_mixed_cell.dual_cells)
        new_mixed_cell.dual_cells[i].dualVector = new_mixed_cell.dual_cells[i].dualVector .* TT.(epsilon * [x.data for x in dual_direction_vectors[i]])
    end

    final_point = stable_intersection_point(new_mixed_cell)

    return final_point - initial_point
end