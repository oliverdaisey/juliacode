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

import Oscar.tropical_pluecker_vector

r"""
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

r"""
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