include("subdivision_of_points.jl")

function regular_subdivisions(s::MixedCell)
    regularSubdivisions = []
    for dualCell in dual_cells(s)
        push!(regularSubdivisions, subdivision_of_points_workaround(points(ambient_support(dualCell)), dual_vector(dualCell)))
    end

    return regularSubdivisions
end

function regular_subdivisions(s::MixedCell, mixedVector::Vector{<:Vector{<:Oscar.TropicalSemiringElem}})
    regularSubdivisions = []
    for (dualCell, dualVector) in zip(dual_cells(s), mixedVector)
        push!(regularSubdivisions, subdivision_of_points_workaround(points(ambient_support(dualCell)), dualVector))
    end

    return regularSubdivisions
end