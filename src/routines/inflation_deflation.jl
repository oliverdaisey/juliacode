
function inflation(dualCell::DualCell, breakingPoint::Vector{Oscar.TropicalSemiringElem{minOrMax}}, drift::Vector{QQFieldElem}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    
    # create all dual cells with respect to the breaking point
    breakingPointDualCells = dual_cells(ambient_support(dualCell), breakingPoint)

    epsilon = QQ(1//1000)

    inflationDualCells = DualCell[]

    # caution: it is possible to get stuck in this loop
    while length(inflationDualCells) == 0
        epsilon *= epsilon # in case epsilon is not sufficiently small
        driftedCenter = center(dualCell) - epsilon*drift
        for s in breakingPointDualCells
            if driftedCenter in convex_hull(s)
                push!(inflationDualCells, s)
            end
        end
    end

    @assert length(inflationDualCells) == 1 "The inflation dual cell is not unique."
    return inflationDualCells[1]

end

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