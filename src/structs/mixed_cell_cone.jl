using Oscar

include("../type_aliases.jl")
"""
Struct encapsulating the defining hyperplanes of a cone.

Each element of `coefficients` is the vector of coefficients ω for the inequality ω⋅x ≤ 0.

Taking the intersection of all such inequalities yields the cone.
"""
struct MixedCellCone

    coefficients::Vector{Vector{QQFieldElem}} # vectors dual to the defining hyperplanes

end

struct CayleyCircuit
    # coefficients defining the corresponding facet
    coefficients::Vector{QQFieldElem}

    # indices of the columns supporting this circuit
    support::Vector{Int}

    # configuqation and column index in cayley embedding giving rise to this circuit
    i::Int
    gamma::Int
end

"""
    mixed_cell_cone(s::GeneralizedMixedCell, M::Matrix{QQFieldElem})

Returns the generalized mixed cell cone of ``s`` in the cayley embedding given by `M`.

INPUTS:
- ``s::GeneralizedMixedCell``: a generalized mixed cell of the cayley embedding. The first component is taken from a loopless facet of the matroid polytope associated to the `PluckerVector`, and components thereafter are from the tropical hypersurfaces.
- ``M::Matrix{QQFieldElem}``: the cayley matrix of the point configurations.
- ``n::Int``: the dimension of the ambient space.
"""
function mixed_cell_cone(s::GeneralizedMixedCell, M::QQMatrix)::MixedCellCone

    # get all indices appearing in generalized mixed cell
    mixed_cell_indices = reduce(vcat, s)
    extra_indices = [i for i in 1:ncols(M) if !(i in mixed_cell_indices)]
    coefficient_vects = Vector{Vector{QQFieldElem}}()

    for extra_index in extra_indices
        push!(coefficient_vects, cone_coefficients(mixed_cell_indices, extra_index, M))
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

    if null_space[nrows(null_space), 1] < 0
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