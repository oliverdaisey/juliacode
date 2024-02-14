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

tropical_pluecker_vector(pluecker_indices::Vector{Vector{Int}}, pluecker_entries::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where minOrMax<:Union{typeof(min),typeof(max)} = TropicalPlueckerVector(pluecker_indices, pluecker_entries)

# TODO: look into printing of tropical pluecker vectors
