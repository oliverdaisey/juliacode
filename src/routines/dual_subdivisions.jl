include("subdivision_of_points.jl")

function dual_subdivisions(s::MixedCell)
    dualSubdivisions = []
    for dualCell in dual_cells(s)
        push!(dualSubdivisions, subdivision_of_points_workaround(ambient_support(dualCell).points, dual_vector(dualCell)))
    end

    return dualSubdivisions
end

function dual_subdivisions(s::MixedCell, mixedVector::Vector{<:Vector{<:Oscar.TropicalSemiringElem}})
    dualSubdivisions = []
    for (dualCell, dualVector) in zip(dual_cells(s), mixedVector)
        push!(dualSubdivisions, subdivision_of_points_workaround(ambient_support(dualCell).points, dualVector))
    end

    return dualSubdivisions
end