using Oscar

"""
Struct encapsulating the defining hyperplanes of a cone.

Each element of `inequalities` is the vector of coefficients ω for the inequality ω⋅x ≤ 0.

Taking the intersection of all such inequalities yields the cone.
"""
struct MixedCellCone

    inequalities::Vector{Vector{QQFieldElem}} # vectors dual to the defining hyperplanes

end

"""
    mixed_cell_cone(s::MixedCell})::MixedCellCone

Returns the mixed cell cone of ``s``.

"""
function mixed_cell_cone(s::MixedCell)::MixedCellCone

    # TODO: Optimise this code, we do not want the cayley embedding explicitly

    # take cayley embedding of m
    M = cayley_embedding([matrix(QQ, points(ambient_support(dual_cells(s)[i]))) for i in 1:length(dual_cells(s))])

    activeIndices = Vector{Int}[]
    offset = 0
    for dualCell in dual_cells(s)
        push!(activeIndices, offset .+ active_indices(dualCell))
        offset += length(active_indices(dualCell))
    end

    # reduce active_indices to a single vector
    activeIndices = vcat(activeIndices...)
    extraIndices = [i for i in 1:nrows(M) if !(i in activeIndices)]
    coefficientVects = Vector{Vector{QQFieldElem}}()

    # {active_indices, extra_indices} form a partition of all indices of the supports
    # each extra index corresponds to (possibly many) hyperplanes defining the cone
    for extraIndex in extraIndices
        push!(coefficientVects, cone_inequalities(activeIndices, extraIndex, transpose(M))...)
    end

    C = MixedCellCone(coefficientVects)

    return C

end

"""
    mixed_cell_cone(mixedCellIndices::Vector{Int}, extra_index::Int, M::QQMatrix)

Returns the inequalities of the definining hyperplane for the mixed cell cone indexed by `extraIndex`.

INPUTS:
- ``mixedCellIndices::Vector{Int}``: the indices of the columns of the cayley matrix defining the generalized mixed cell.
- ``extraIndex::Int``: an extra index not belonging to the mixed cell.
- ``M::QQMatrix``: the cayley matrix of the point configurations, encoded as columns.
"""
function cone_inequalities(mixedCellIndices::Vector{Int}, extraIndex::Int, M::QQMatrix)
    # get submatrix defined by indices in I
    I = copy(mixedCellIndices)
    push!(I, extraIndex)
    submatrix = M[:,I]
    # compute the null space
    ν, nullSpace = Oscar.nullspace(submatrix)
    # @assert ν == 1 "The null space should be one dimensional" # This is actually only true for the hypersurface case
    println("Null space: ", nullSpace)
    coefficients = [nullSpace[:,i] for i in 1:ncols(submatrix)]
    println("Coefficients: ", coefficients)
    
    if nullSpace[ncols(nullSpace), 1] > 0
        # we want the normal to point in the right direction
        nullSpace = -nullSpace
    end

    coefficientDict = Dict{Int, QQFieldElem}(
        i => 0 for i in 1:ncols(M)
    )

    for i in 1:length(I)
        coefficientDict[I[i]] = nullSpace[i]
    end

    return [coefficientDict[i] for i in 1:ncols(M)]
end

function mixed_cell_cone_to_polyhedron(C::MixedCellCone)::Polyhedron

    coefficient_matrix = matrix(QQ, vcat(C.inequalities))
    
    return polyhedron(coefficient_matrix, zeros(QQ, nrows(coefficient_matrix)))
    
end

function inequalities(C::MixedCellCone)::Vector{Vector{QQFieldElem}}
    return C.inequalities
end

function Base.show(io::IO, C::MixedCellCone)
    print(io, "Mixed cell cone with facet inequalities $(inequalities(C))")
end