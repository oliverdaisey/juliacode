"""
    tropical_homotopy_continuation(T::MixedCellTracker)

Runs the tropical homotopy continuation algorithm on the mixed cell tracker `T` and returns all the mixed cells at the target dual weight. This works as a depth-first search on the tree of all mixed cell trackers.
"""
function tropical_homotopy_continuation(T::MixedCellTracker)
    # begin with a deep copy of T to ensure no side effects on the original tracker
    T = deepcopy(T)
    mixedCellsAtTarget = MixedCell[]
    find_mixed_cells(T, mixedCellsAtTarget)
    return mixedCellsAtTarget
end

function find_mixed_cells(T::MixedCellTracker, mixedCellsAtTarget::Vector{MixedCell} = MixedCell[])

    # if we have a cell to visit, proceed
    if length(pointers(T)) > 1
        pointOfInterest, changingSupports = next_point_of_interest(T)
        # if no supports change, then move to the next node
        if changingSupports == []
            return find_mixed_cells(T, mixedCellsAtTarget)
        end

        # otherwise we are in the case that supports change, so we need to calculate new trackers
        tropicalDrift = tropical_drift(T)

        deflationsVector = Vector{DualCell}[]
        for j in 1:length(dual_cells(mixed_cell(T)))
            dualCell = dual_cells(mixed_cell(T))[j]
            inflatedDualCell = inflation(dualCell, pointOfInterest[j], tropicalDrift)
            deflatedVector = deflation(inflatedDualCell, tropicalDrift)

            # if we have a deflation, this support is changing
            if length(deflatedVector) > 0
                push!(deflationsVector, deflatedVector)
            else
                # otherwise the support does not change
                push!(deflationsVector, [dualCell])
            end
        end
        
        # now we have all the new dual cells, combine them into all possible mixed cells at this point, and get the new trackers
        # each mixed cell will use exactly one dual cell from each deflationVector in deflationsVector

        # calculate the number of indices
        numIndices = length(deflationsVector)

        # create an array of indices to iterate over
        indices = [1:length(deflationsVector[i]) for i in 1:numIndices]

        # generate all possible combinations of indices
        combinations = Iterators.product(indices...)

        # create an empty array to store the new mixed cells
        newMixedCells = []

        # iterate over the combinations
        for combination in combinations
            # create an empty array to store the dual cells for this combination
            dualCells = DualCell[]
            
            # iterate over the indices in the combination
            for i in 1:numIndices
                # get the dual cell at the corresponding index
                dualCell = deflationsVector[i][combination[i]]
                
                # add the dual cell to the array
                push!(dualCells, dualCell)
            end
            
            # create a new mixed cell using the dual cells
            mixedCell = mixed_cell(dualCells)
            
            # add the new mixed cell to the array
            push!(newMixedCells, mixedCell)
        end

        
        # now generate all new trackers from the new mixed cells
        newTrackers = [mixed_cell_tracker(deepcopy(mixed_path(T)), deepcopy(mixedCell)) for mixedCell in newMixedCells]

        # run the continuation on each new tracker
        for newTracker in newTrackers
            find_mixed_cells(newTracker, mixedCellsAtTarget)
        end
        

    else
        # otherwise this tracker has reached the end of its course, do something...
        push!(mixedCellsAtTarget, mixed_cell(T))
    end
end