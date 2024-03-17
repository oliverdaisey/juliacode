using Oscar

"""
Struct encapsulating the defining hyperplanes of a cone.

Each element of `coefficients` is the vector of coefficients ω for the inequality ω⋅x ≤ 0.

Taking the intersection of all such inequalities yields the cone.
"""
struct MixedCellCone

    coefficients::Vector{Vector{QQFieldElem}} # vectors dual to the defining hyperplanes

end

"""
    mixed_cell_cone(s::GeneralizedMixedCell, M::Matrix{QQFieldElem})

Returns the mixed cell cone of ``s``.

"""
function mixed_cell_cone(s::MixedCell)::MixedCellCone

    # TODO: Optimise this code, we do not want the cayley embedding explicitly

    # take cayley embedding of m
    M = cayley_embedding([matrix(QQ, s.dual_cells[i].ambientDualSupport.points) for i in 1:length(s.dual_cells)])

    activeIndices = Vector{Int}[]
    offset = 0
    for dual_cell in s.dual_cells
        push!(activeIndices, offset .+ active_indices(dual_cell))
        offset += length(active_indices(dual_cell))
    end

    # reduce active_indices to a single vector
    activeIndices = vcat(activeIndices...)
    extraIndices = [i for i in 1:nrows(M) if !(i in activeIndices)]
    coefficient_vects = Vector{Vector{QQFieldElem}}()

    # {active_indices, extra_indices} form a partition of all indices of the supports
    # each extra index corresponds to a facet of the mixed cell cone
    for extraIndex in extraIndices
        push!(coefficient_vects, cone_coefficients(activeIndices, extraIndex, transpose(M))) # transpose to fit previous implementation
    end

    C = MixedCellCone(coefficient_vects)

    return C

end

"""
    mixed_cell_cone(mixed_cell_indices::Vector{Int}, extra_index::Int, M::QQMatrix)

Returns the coefficients of the definining hyperplane for the mixed cell cone indexed by `extra_index`.

INPUTS:
- ``mixed_cell_indices::Vector{Int}``: the indices of the columns of the cayley matrix defining the generalized mixed cell.
- ``extra_index::Int``: an extra index not belonging to the mixed cell.
- ``M::QQMatrix``: the cayley matrix of the point configurations, encoded as columns.
"""
function cone_coefficients(mixed_cell_indices::Vector{Int}, extra_index::Int, M::QQMatrix)
    # get submatrix defined by indices in I
    I = copy(mixed_cell_indices)
    push!(I, extra_index)
    submatrix = M[:,I]
    # compute the null space
    ν, null_space = Oscar.nullspace(submatrix)
    @assert ν == 1 "The null space should be one dimensional"

    if null_space[ncols(null_space), 1] > 0
        # we want the normal to point in the right direction
        null_space = -null_space
    end

    coefficient_dict = Dict{Int, QQFieldElem}(
        i => 0 for i in 1:ncols(M)
    )

    for i in 1:length(I)
        coefficient_dict[I[i]] = null_space[i]
    end

    return [coefficient_dict[i] for i in 1:ncols(M)]
end

function mixed_cell_cone_to_polyhedron(C::MixedCellCone)::Polyhedron

    coefficient_matrix = matrix(QQ, vcat(C.coefficients))
    
    return polyhedron(coefficient_matrix, zeros(QQ, nrows(coefficient_matrix)))
    
end