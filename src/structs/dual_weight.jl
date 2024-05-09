MinOrMax = Union{typeof(min),typeof(max)}

"""
The dual type is important for computing perturbations.
"""
mutable struct DualWeight{dualType<:DualType, minOrMax<:Union{typeof(min), typeof(max)}}

    entries::Dict{Vector{Int}, Oscar.TropicalSemiringElem{minOrMax}} 
    data # to do: add type for this, it's either a polynomial, matrix + valuation, or empty

    function DualWeight{dualType, minOrMax}(data) where {dualType<:DualType, minOrMax<:Union{typeof(min), typeof(max)}}
        return new{dualType, minOrMax}(Dict{Vector{Int}, Oscar.TropicalSemiringElem{minOrMax}}(), data)
    end
end

function dual_type(w::DualWeight)
    return typeof(w).parameters[1]
end

function cache(w::DualWeight, index::Vector{Int}, weight::Oscar.TropicalSemiringElem{minOrMax}) where minOrMax<:MinOrMax
    w.entries[index] = weight
end

function data(w::DualWeight)
    return w.data
end 

function Base.show(io::IO, w::DualWeight)
    print(io, "Dual weight of type $(dual_type(w))")
end

"""
    dual_weight(polynomial::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{minOrMax}}) where minOrMax<:Union{typeof(min), typeof(max)}

    Create a Hypersurface dual weight from a polynomial.
"""
function dual_weight(polynomial::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{minOrMax}}) where {minOrMax<:Union{typeof(min), typeof(max)}}
    return DualWeight{Hypersurface, minOrMax}(polynomial)
end

function dual_weight(dualType::LinearType, realisation::MatElem{K}, nu::TropicalSemiringMap{L,t,minOrMax}) where {K,L,t,minOrMax,LinearType<:Union{typeof(Linear), typeof(InvertedLinear)}}
    return DualWeight{dualType, minOrMax}((realisation, nu))
end

"""
    dual_weight(dualType::LinearType, plueckerIndices::Vector{Vector{Int}}, plueckerEntries::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where {minOrMax<:Union{typeof(min), typeof(max)}, LinearType<:Union{typeof(Linear), typeof(InvertedLinear)}}

Create a LinearType dual weight from a predefined active support and weights. This is useful in the case that your dual weight is not realisable in small dimension (or you just don't have a realisation to hand.)

"""
function dual_weight(dualType::LinearType, activeSupport::Vector{Vector{Int}}, weights::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where {minOrMax<:Union{typeof(min), typeof(max)}, LinearType<:Union{typeof(Linear), typeof(InvertedLinear)}}

    @assert length(activeSupport) == length(weights) "The number of indices and entries must match"
    @assert dualType == Linear || all(issubset(index, [0,1]) for index in activeSupport) "Active support must be a 0/1 vectors for an InvertedLinear dual weight"
    @assert dualType == InvertedLinear || all(issubset(index, [0,-1]) for index in activeSupport) "Active support must be a 0/-1 vectors for a Linear dual weight"

    # construct a dual weight with no data
    w = DualWeight{dualType, minOrMax}(nothing)
    # initialise the given entries
    for i in 1:length(activeSupport)
        cache(w, activeSupport[i], weights[i])
    end

    return w
end

function Base.getindex(w::DualWeight, index::Vector{Int})
    if haskey(w.entries, index)
        return w.entries[index]
    end
    return compute_entry(w, index)
end

function compute_entry(w::DualWeight{Hypersurface, minOrMax}, index::Vector{Int}) where minOrMax<:MinOrMax
    weight = coeff(data(w), index)
    cache(w, index, weight)
    return weight
end

function compute_entry(w::DualWeight{LinearType, minOrMax}, index::Vector{Int}) where {minOrMax<:MinOrMax, LinearType<:Union{Linear, InvertedLinear}}

    # assert that data exists
    @assert !isempty(data(w)) "Data must be set for this dual weight."

    realisation, nu = data(w)

    # perform some checks
    @assert length(index) == ncols(realisation) "Index must have the same length as the number of columns in the realisation (expected $(ncols(realisation)), got $(length(index)))"
    @assert dual_type(w) == Linear || issubset(index, [0,1]) "Index must be a 0/1 vector for an InvertedLinear dual weight"
    @assert dual_type(w) == InvertedLinear || issubset(index, [0,-1]) "Index must be a 0/-1 vector for a Linear dual weight"
    

    weight = nu(det(realisation[:,findall(!iszero, index)]))
    cache(w, index, weight)
    return weight
end

function Base.:/(w1::DualWeight, w2::DualWeight)
    error("Division of dual weights is not currently supported.")
end