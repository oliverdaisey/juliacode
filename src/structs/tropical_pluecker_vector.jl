using Oscar
import LinearAlgebra.rank

struct TropicalPlueckerVector{minOrMax<:Union{typeof(min),typeof(max)}}
    pluecker_indices::Vector{Vector{Int}}
    pluecker_entries::Vector{Oscar.TropicalSemiringElem{minOrMax}}

    function TropicalPlueckerVector(pluecker_indices::Vector{Vector{Int}}, pluecker_entries::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where minOrMax<:Union{typeof(min),typeof(max)}
        @assert length(pluecker_indices)==length(pluecker_entries) "The number of indices and entries must match"
        @assert minOrMax==typeof(min) "max tropical pluecker vectors currently unsupported"
        return new{minOrMax}(pluecker_indices, pluecker_entries)
    end
end

struct LazyTropicalPlueckerVector{minOrMax<:Union{typeof(min), typeof(max)}}
    realisation::Union{Matrix{QQFieldElem}, QQMatrix}
    nu::TropicalSemiringMap{Kt,t,minOrMax} where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem}
    tropicalPlueckerVector::Dict{Vector{Int}, Oscar.TropicalSemiringElem{minOrMax}}

    function LazyTropicalPlueckerVector(realisation::Union{Matrix{QQFieldElem}, QQMatrix}, nu::TropicalSemiringMap{Kt,t,minOrMax}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
        return new{minOrMax}(realisation, nu, Dict{Vector{Int}, Oscar.TropicalSemiringElem{minOrMax}}())
    end
end

function rank(T::LazyTropicalPlueckerVector)
    return size(T.realisation, 1)
end

function Base.getindex(T::LazyTropicalPlueckerVector, index::Vector{Int})
    if haskey(T.tropicalPlueckerVector, index)
        return T.tropicalPlueckerVector[index]
    else
        # check the input was valid
        @assert length(index)==rank(T) "The index must have the same length as the rank of the tropical pluecker vector (expected $(rank(T)), got $(length(index)))"
        # compute the entry
        entry = nu(det(T.realisation[:, index]))
        T.tropicalPlueckerVector[index] = entry
        return entry
    end
end

import Oscar.tropical_pluecker_vector

"""
    tropical_pluecker_vector(pluecker_indices::Vector{Vector{Int}}, pluecker_entries::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where minOrMax<:Union{typeof(min),typeof(max)}

Create a tropical Pluecker vector from the given indices and entries.

# Arguments
- `pluecker_indices::Vector{Vector{Int}}`: The indices of the pluecker vector.
- `pluecker_entries::Vector{Oscar.TropicalSemiringElem{minOrMax}}`: The entries of the pluecker vector. These are ordered according to the indices.
- `minOrMax::Union{typeof(min),typeof(max)}`: The min or max convention.

# Returns
A tropical pluecker vector from the given indices and entries.

"""
tropical_pluecker_vector(pluecker_indices::Vector{Vector{Int}}, pluecker_entries::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where minOrMax<:Union{typeof(min),typeof(max)} = TropicalPlueckerVector(pluecker_indices, pluecker_entries)

function Base.show(io::IO, T::TropicalPlueckerVector)
    # make sure we can print entries withotu type information
    entries_string = ""
    for i in range(1, (length(T.pluecker_entries)-1))
        entries_string *= string(T.pluecker_entries[i]) * ", "
    end
    entries_string *= string(T.pluecker_entries[end])

    print(io, "Tropical Pluecker vector on indices $(T.pluecker_indices) with entries [$entries_string]")
end

"""
    rank(T::TropicalPlueckerVector)

Compute the rank of a tropical pluecker vector. This is the rank of the matroid with bases given by the indices.

# Arguments
- `T::TropicalPlueckerVector`: The tropical pluecker vector.

# Returns
The rank of the tropical pluecker vector.

"""
function rank(T::TropicalPlueckerVector)
    # this is the number of unique integers appearing in the pluecker indices
    return length(unique(vcat(T.pluecker_indices...)))
end

function pluecker_indices(T::TropicalPlueckerVector)
    return T.pluecker_indices
end

function pluecker_entries(T::TropicalPlueckerVector)
    return T.pluecker_entries
end

function matroid_polytope_vertices(n::Int, k::Int)
    subsets = Oscar.subsets([i for i in 1:n], k)
    indicator_matrix = zeros(Int, length(subsets), n)
    for (i, subset) in enumerate(subsets)
        for index in subset
            indicator_matrix[i, index] = 1
        end
    end
    return indicator_matrix
end