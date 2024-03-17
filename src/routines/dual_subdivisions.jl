
function regular_subdivisions(s::MixedCell)
    regularSubdivisions = []
    for dualCell in dual_cells(s)
        push!(regularSubdivisions, subdivision_of_points_workaround(points(ambient_support(dualCell)), dual_vector(dualCell)))
    end

    return regularSubdivisions
end

function regular_subdivisions(s::MixedCell, mixedVector::Vector{<:Vector{<:Oscar.TropicalSemiringElem}})
    @assert length(dual_cells(s)) == length(mixedVector) "The number of dual cells and the number of mixed vectors do not match."
    regularSubdivisions = []
    for (dualCell, dualVector) in zip(dual_cells(s), mixedVector)
        @assert length(dualVector) == length(ambient_support(dualCell)) "The length of the dual vector does not match the number of points in the dual cell."
        push!(regularSubdivisions, subdivision_of_points_workaround(points(ambient_support(dualCell)), dualVector))
    end

    return regularSubdivisions
end