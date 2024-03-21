
function tropical_homotopy_continuation(T::MixedCellTracker)
    # if we have a cell to visit, proceed
    if length(pointers(T)) > 1
        pointOfInterest, changingSupports = next_point_of_interest(T)

        # if no supports change, then move to the next node
        # remove the head of the pointers, and call this with the new mixed cell tracker
        tropicalDrift = tropical_drift(T)

        deflationsVector = Vector{DualCell}[]
        for j in 1:length(dual_cells(mixed_cell(T)))
            dualCell = dual_cells(mixed_cell(T))[j]
            inflatedDualCell = inflation(dualCell, pointOfInterest[j], tropicalDrift)
            deflatedVector = deflation(inflatedDualCell, drift)

            # if we have a deflation, this support is changing
            if length(deflatedVector) > 0
                push!(deflationsVector, deflatedVector)
            else
                # otherwise the support does not change
                push!(deflationsVector, [dualCell])
            end
        end
        
        # now we have all the new dual cells, combine them into all possible mixed cells at this point, and get the new trackers
        # TODO: write this. Then with all the new mixed cell trackers, run this function recursively
    else
        # otherwise this tracker has reached the end of its course, do something...
        println(mixed_cell(T))
        # TODO: make this do something more useful (yield?)
    end
end