using Oscar

"""
ω is generic means it is not in the codimension 2 skeleton of the secondary fan of
cayley_embedding
"""

"""
Type definitions
"""
MixedCell = Vector{NTuple{2,Int}}
PointConfiguration = Vector{Matrix{QQFieldElem}}
GeneralizedMixedCell = Vector{Vector{Int}}
# alias for a tuple of tropical polynmomials
TropicalTuple{N} = NTuple{N, AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}}
TropicalPoly = AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}} # alias for a tropical polynomial
PluckerVector = Pair{Vector{Vector{Int}}, Vector{Int}}

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
End type definitions
"""


"""
Utility to remove an item from a vector
"""
function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

"""
Function to compute the bistellar flip as described in Anders' paper.

INPUTS:
- `A::Vector{PointConfiguration}`: a vector of point configurations in the same dimension
- `σ::MixedCell`: a mixed cell of the cayley embedding. We encode it as a choice of two indices from each
column of the cayley embedding in each dimension.
- `c::Vector{QQFieldElem}`: a circuit defining a facet F of the mixed cell cone C_σ
- `A_mix::Vector{MixedCell}`: a vector of mixed cells sufficiently close to F
OUTPUTS:
- `B_mix::Vector{MixedCell}`: a vector of mixed cells sufficiently close to F after the bistellar flip

"""
function bistellar_flip(
    A::Vector{PointConfiguration}, σ::MixedCell, c::CayleyCircuit, A_mix::Vector{MixedCell}
)::Vector{MixedCell}

    # step 1
    B_mix = Vector{MixedCell}()

    # step 2 is handled by the caller
    i = c.I
    gamma = c.gamma

    # step 3
    α, β = σ[i]

    # step 4
    for τ in A_mix
        if gives_rise(τ, c)
            mixed_cell_cayley_indices = reduce(vcat, σ)
            if c.coefficients[α] > 0
                updated_indices = deepcopy(mixed_cell_cayley_indices)
                push!(updated_indices, gamma)
                remove!(updated_indices, α)
                if updated_indices ∉ B_mix
                    push!(B_mix, updated_indices)
                end 
            end
            # note we can update B_mix twice in one step!
            if c.coefficients[β] > 0
                updated_indices = deepcopy(mixed_cell_cayley_indices)
                push!(updated_indices, gamma)
                remove!(updated_indices, β)
                if updated_indices ∉ B_mix
                    push!(B_mix, updated_indices)
                end 
            end
        end
    end

    # step 5
    return B_mix
    
end

"""
Decides whether a mixed cell τ gives rise to a circuit
"""
function gives_rise(τ::MixedCell, c::CayleyCircuit)::Bool
    
    
    
end


"""
    ray_intersects_cone(C::Cone, u::Vector{QQFieldElem}, v::Vector{QQFieldElem})

Determines the value of ``t`` at which the ray ``u + tv`` with ``t in [0, 1]`` first intersects the boundary of the cone ``C``.

In our problem we will have ``v = w_target - w_start`` as we try to walk the shortest path to the target.

INPUTS:
- ``C::Cone``: cone
- ``u::Vector{QQFieldElem}``: start point of the ray, required to be in C
- ``v::Vector{QQFieldElem}``: direction of the ray, required to be nonzero

OUTPUTS:
- ``t::QQFieldElem``: parameter value at which the ray first intersects the boundary of the cone, or ``nothing`` if the ray does not intersect the cone.

# Examples
```jldoctest
julia> C = cone_from_inequalities([-1 0 0; 0 -1 0; 0 0 -1])
julia> u = Vector{QQFieldElem}([1, 1, 1])
julia> v = Vector{QQFieldElem}([-1, 0, 0])
julia> ray_intersects_cone(C, u, v)
1
julia> v[1] = 1
julia> ray_intersects_cone(C, u, v)
nothing
```
"""
function ray_intersects_cone(C::Cone, u::Vector{QQFieldElem}, v::Vector{QQFieldElem})

    @req u in C "u must be in the cone C"
    @req !iszero(v) "v must be nonzero"
    
    intersection_points = []
    for F in facets(C)
        v_dot = dot(F.a , v)
        if v_dot == 0
            # this means they are parallel, so we can't intersect
            continue
        end
        t = dot(-F.a , u) / v_dot
        if t >= 0
            push!(intersection_points, t)
        end
    end
    
    if isempty(intersection_points)
        # nothing satisfied the criteria, so ray does not intersect
        return nothing
    else
        # intersection point is precisely the minimum of the intersection points
        return minimum(intersection_points)
    end


end

"""
    mixed_cell_cone(σ::MixedCell, cayley_matrix::QQMatrix)

Returns the mixed cell cone of ``σ`` in the cayley embedding given by `cayley_matrix`.

INPUTS:
- ``σ::MixedCell``: a mixed cell of the cayley embedding. We encode it as a choice of two indices from each
column of the cayley embedding in each dimension.
- ``cayley_matrix::QQMatrix``: the cayley matrix of the point configuration

OUTPUTS:
- ``C_σ::Cone``: the mixed cell cone of ``σ`` in the cayley embedding given by `cayley_matrix`.
"""
function mixed_cell_cone(σ::MixedCell, cayley_matrix::QQMatrix)

    # we implement this by taking the facet normals from the null space of the cayley matrix
    
    # first expand σ into a vector of ints
    σ_vector = Vector{Int}()
    for tuple in σ
        push!(σ_vector, tuple...)
    end

    # each element of this array corresponds to a facet
    indices_of_cayley = [i for i in 1:nrows(cayley_matrix) if !(i in σ_vector)]
    circuits = []

    for γ in indices_of_cayley

        # ν nullity, null_space basis for null null_space
        # TODO: make this work correctly
        ν, null_space = Oscar.nullspace(cayley_matrix[reduce(vcat, [σ_vector,γ]),:])

        println(null_space)
        first_vector = findfirst(x->x[γ] != 0, eachcol(Matrix(null_space)))

        if first_vector[γ] > 0
            # we want the normal to point inwards
            first_vector = -first_vector
        end

        # TODO: need to pad out first_vector with zeroes according to the other indices

        # this is our circuit for the mixed cell cone
        push!(circuits, first_vector)
    end

    # at this point we use the circuits to prduce the mixed cell cone
    return cone_from_inequalities(Matrix(transpose(hcat(circuits...))))

    
end

"""
    mixed_cell_cone(s::GeneralizedMixedCell, M::Matrix{QQFieldElem})

Returns the generalized mixed cell cone of ``s`` in the cayley embedding given by `M`.

INPUTS:
- ``s::GeneralizedMixedCell``: a generalized mixed cell of the cayley embedding. The first component is taken from a loopless facet of the matroid polytope associated to the `PluckerVector`, and components thereafter are from the tropical hypersurfaces.
- ``M::Matrix{QQFieldElem}``: the cayley matrix of the point configurations.
- ``n::Int``: the dimension of the ambient space.
"""
function mixed_cell_cone(s::GeneralizedMixedCell, M::QQMatrix)

    # get all indices appearing in generalized mixed cell
    mixed_cell_indices = reduce(vcat, s)
    extra_indices = [i for i in 1:ncols(M) if !(i in mixed_cell_indices)]
    coefficient_vects = Vector{Vector{QQFieldElem}}()

    for extra_index in extra_indices
        push!(coefficient_vects, cone_coefficients(mixed_cell_indices, extra_index, M))
    end

    coefficient_matrix = matrix(QQ, vcat(coefficient_vects))
    
    return polyhedron(coefficient_matrix, zeros(QQ, nrows(coefficient_matrix)))

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