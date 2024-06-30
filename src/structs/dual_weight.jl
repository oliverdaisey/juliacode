MinOrMax = Union{typeof(min),typeof(max)}

"""
    DualWeight{dualType<:DualType, minOrMax<:Union{typeof(min), typeof(max)}}

A dual weight is a choice of coefficients/weights for a dual cell. From dual weights one can construct dual paths, mixed paths, and compute stable intersections. Additionally parametrises a dual type, which provides constraints for computing perturbations and deflations. 

# Fields
- `entries::Dict{Vector{Int}, Oscar.TropicalSemiringElem{minOrMax}}`: A dictionary of indices to weights. Allows for caching of computed weights.
- `data`: The data from which entries of the dual weight may be computed. This is either a polynomial (`Hypersurface`), a matrix and valuation (`Linear`/`InvertedLinear`), or nothing.
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

"""
    dual_weight(dualType::LinearType, realisation::MatElem{K}, nu::TropicalSemiringMap{L,t,minOrMax}) where {K,L,t,minOrMax,LinearType<:Union{typeof(Linear), typeof(InvertedLinear)}}

Create a LinearType dual weight from a realisation matrix and a tropical semiring map. Entries of the dual weight are computed on the fly and cached.
"""
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

"""
    compute_entry(w::DualWeight{LinearType, minOrMax}, index::Vector{Int}) where {minOrMax<:MinOrMax, LinearType<:Union{Linear, InvertedLinear}}

We index linear dual weights by their indicator vectors (i.e. vertices of the corresponding matroid polytope).
"""
function compute_entry(w::DualWeight{LinearType, minOrMax}, index::Vector{Int}) where {minOrMax<:MinOrMax, LinearType<:Union{Linear, InvertedLinear}}

    if isnothing(data(w))
        println("w has entries $(w.entries)")
        return w.entries[index]
    end

    realisation, nu = data(w)

    # perform some checks
    @assert length(index) == ncols(realisation) "Index must have the same length as the number of columns in the realisation (expected $(ncols(realisation)), got $(length(index)))"
    @assert dual_type(w) == Linear || issubset(index, [0,1]) "Index must be a 0/1 vector for an InvertedLinear dual weight"
    @assert dual_type(w) == InvertedLinear || issubset(index, [0,-1]) "Index must be a 0/-1 vector for a Linear dual weight"
    # make sure the correct number of entries is nonzero

    weight = nu(det(realisation[:,findall(!iszero, index)]))
    cache(w, index, weight)
    return weight
end

function Base.length(w::DualWeight{Hypersurface, minOrMax}) where minOrMax<:MinOrMax
    return length(w.entries)
end

function Base.iterate(w::DualWeight{Hypersurface, minOrMax}) where minOrMax<:MinOrMax
    return iterate(w.entries)
end

function Base.iterate(w::DualWeight{Hypersurface, minOrMax}, state) where minOrMax<:MinOrMax
    return iterate(w.entries, state)
end

"""
    Base.:/(w1::DualWeight, w2::DualWeight)

Division of dual weights is not supported, and likely does not make sense.
"""
function Base.:/(w1::DualWeight, w2::DualWeight)
    error("Division of dual weights is not currently supported.")
end

"""
    direction(w1::DualWeight{Hypersurface, minOrMax}, w2::DualWeight{Hypersurface, minOrMax})

Compute the direction between two Hypersurface dual weights.
"""
function direction(w1::DualWeight{Hypersurface, minOrMax}, w2::DualWeight{Hypersurface, minOrMax}) where minOrMax<:MinOrMax

    i1 = collect(indices(w1))
    i2 = collect(indices(w2))
    @assert i1 == i2 "The ambient dual support changed between the two dual weights"

    # compute the direction by evaluating both dual weights entirely, then taking the difference
    entries = []
    for index in i1
        push!(entries, getindex(w1, index) / getindex(w2, index))
    end

    return entries

end

"""
    direction(w1::DualWeight{LinearType, minOrMax}, w2::DualWeight{LinearType, minOrMax})

Compute the direction between two LinearType dual weights.
"""
function direction(w1::DualWeight{LinearType, minOrMax}, w2::DualWeight{LinearType, minOrMax}) where {minOrMax<:MinOrMax, LinearType<:Union{Linear, InvertedLinear}}

    error("Computation of direction for linear dual weights is not yet implemented")
    println("Warning: slow computation of direction for linear dual weights")
    # compute the direction by evaluating both dual weights entirely, then taking the difference
    entries = []

    
end

function indices(w::DualWeight{Hypersurface, minOrMax}) where minOrMax<:MinOrMax
    
    # get all exponent vectors from data
    return exponents(data(w))
end

function indices(w::DualWeight{LinearType, minOrMax}) where {minOrMax<:MinOrMax, LinearType<:Union{Linear, InvertedLinear}}

    # if we do not have any data, it means that all the entries are cached
    if isnothing(data(w))
        return sort(collect(keys(w.entries)))
    end

    # otherwise we must have a realisation
    n = nrows(data(w))
    m = ncols(data(w))
    
    return Oscar.combinations(n, m)
end

function Base.length(w::DualWeight{Hypersurface, minOrMax}) where minOrMax<:MinOrMax
    return length(indices(w))
end

function Base.length(w::DualWeight{LinearType, minOrMax}) where {minOrMax<:MinOrMax, LinearType<:Union{Linear, InvertedLinear}}

    # if we do not have any data, it means that all the entries are cached
    println("data is $(data(w))")
    if isnothing(data(w))
        return length(keys(w.entries))
    end

    # otherwise we must have a realisation
    n = nrows(data(w))
    m = ncols(data(w))
    
    return Oscar.combinations(n, m)
end

import Base.+
function (+)(w1::DualWeight{Hypersurface, minOrMax}, v::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where minOrMax<:MinOrMax

    @assert length(v) == length(indices(w1)) "Dual weight and vector of tropical numbers is not the same length"

    # create a new dual weight
    w = DualWeight{Hypersurface, minOrMax}(data(w1))
    for (index, weight) in zip(indices(w1), v)
        w.entries[index] = w[index] * weight # tropical multiplciation is addition
    end

    return w
    
end

function +(w1::DualWeight{LinearType, minOrMax}, v::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where {minOrMax<:MinOrMax, LinearType<:Union{Linear, InvertedLinear}}

    @assert length(v) == length(w1) "Dual weight and vector of tropical numbers is not the same length"

    # create a new dual weight
    w = deepcopy(w1)
    for (index, weight) in zip(indices(w1), v)
        if isone(weight)
            w.entries[index] = w[index] * weight # tropical multiplciation is addition
        end
    end

    return w
    
end