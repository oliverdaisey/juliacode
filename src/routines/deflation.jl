
function deflation(inflatedCell::DualCell, drift::Vector{QQFieldElem})
    
    deflatedCells = DualCell[]
    for face in dual_facets(inflatedCell)
        epsilon = QQ(1//1000000)
        isPlusIn = false
        isMinusIn = false
        while !isPlusIn && !isMinusIn
            driftedCenterPlus = center(face) + epsilon*drift
            driftedCenterMinus = center(face) - epsilon*drift
            isPlusIn = (driftedCenterPlus in convex_hull(inflatedCell))
            isMinusIn = (driftedCenterMinus in convex_hull(inflatedCell))
            if isPlusIn
                push!(deflatedCells, face)
            end
            epsilon *= epsilon # in case epsilon is not sufficiently small
        end
    end

    return deflatedCells
end