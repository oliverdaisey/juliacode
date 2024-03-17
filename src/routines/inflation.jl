
function inflation(dualCell::DualCell, breakingPoint::Vector{Oscar.TropicalSemiringElem{minOrMax}}, drift::Vector{QQFieldElem}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    
    # create all dual cells with respect to the breaking point
    breakingPointDualCells = dual_cells(ambient_support(dualCell), breakingPoint)

    epsilon = QQ(1//1000000)
    driftedCenter = center(dualCell) - epsilon*drift

    inflationDualCells = DualCell[]

    while length(inflationDualCells) == 0
        for s in breakingPointDualCells
            if driftedCenter in convex_hull(s)
                push!(inflationDualCells, s)
            end
        end
        epsilon *= epsilon # in case epsilon is not sufficiently small
    end

    @assert length(inflationDualCells) == 1 "The inflation dual cell is not unique."
    return inflationDualCells[1]

end