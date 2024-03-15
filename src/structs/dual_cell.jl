using Oscar
include("support.jl")
include("dual_type.jl")

mutable struct DualCell{cellType<:DualType,minOrMax<:Union{typeof(min),typeof(max)}}

    ambientSupport::Support{cellType}
    activeSupport::Vector{Int}
    dualVector::Vector{Oscar.TropicalSemiringElem{minOrMax}}

    function DualCell{cellType,minOrMax}(ambientSupport::Support{cellType}, activeSupport::Vector{Int}, dualVector::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where {cellType<:DualType,minOrMax<:Union{typeof(min),typeof(max)}}
        return new{cellType,minOrMax}(ambientSupport, activeSupport, dualVector)
    end
end

"""
    dual_cell(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, dualVector::Vector{Oscar.TropicalSemiringElem{minOrMax}}, cellType::Symbol, ::typeof(min)=min)

Create a dual cell of the given type, using the min convention, with given ambient and active support.

# Arguments
- `ambientSupport::Matrix{Int}`: The ambient support of the dual cell, with columns as points.
- `activeSupport::Vector{Int}`: The active support of the dual cell, corresponding to indices of columns of the ambient support.
- `dualVector::Vector{Oscar.TropicalSemiringElem{minOrMax}}`: The dual vector giving rise to this dual cell.
- `cellType::Symbol`: The type of the dual cell, must be one of :hypersurface, :linear, or :inverted_linear.
- `::typeof(min)`: The min convention.

# Returns
A dual cell of the given type, using the min convention, with given ambient and active support.
"""
function dual_cell(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, dualVector::Vector{Oscar.TropicalSemiringElem{typeof(min)}}, cellType::Symbol, ::typeof(min)=min)
    check_dual_cell_inputs(ambientSupport, activeSupport, cellType)
    if cellType == :hypersurface
        return DualCell{Hypersurface,typeof(min)}(ambientSupport, activeSupport)
    elseif cellType == :linear
        return DualCell{Linear,typeof(min)}(ambientSupport, activeSupport)
    elseif cellType == :inverted_linear
        return DualCell{InvertedLinear,typeof(min)}(ambientSupport, activeSupport)
    else
        throw(ArgumentError("cellType must be one of :hypersurface, :linear, or :inverted_linear"))
    end
end

"""
    dual_cell(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, dualVector::Vector{Oscar.TropicalSemiringElem{minOrMax}}, cellType::Symbol, ::typeof(max)=max)

Create a dual cell of the given type, using the max convention, with given ambient and active support.

# Arguments
- `ambientSupport::Matrix{Int}`: The ambient support of the dual cell, with columns as points.
- `activeSupport::Vector{Int}`: The active support of the dual cell, corresponding to indices of columns of the ambient support.
- `dualVector::Vector{Oscar.TropicalSemiringElem{minOrMax}}`: The dual vector giving rise to this dual cell.
- `cellType::Symbol`: The type of the dual cell, must be one of :hypersurface, :linear, or :inverted_linear.
- `::typeof(max)`: The max convention.

# Throws
An error, as max tropical dual cells are currently unsupported.

"""
function dual_cell(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, dualVector::Vector{Oscar.TropicalSemiringElem{typeof(max)}}, cellType::Symbol, ::typeof(max))
    error("max tropical dual cells currently unsupported")
end

"""
    check_dual_cell_inputs(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, cellType::Symbol, pluecker_vector::Union{TropicalPlueckerVector, Nothing}=nothing)

Check that the inputs to create a dual cell are valid.

# Arguments
- `ambientSupport::Matrix{Int}`: The ambient support of the dual cell.
- `activeSupport::Vector{Int}`: The active support of the dual cell, corresponding to indices of columns of the ambient support.
- `cellType::Symbol`: The type of the dual cell, must be one of :hypersurface, :linear, or :inverted_linear.
- `pluecker_vector::Union{TropicalPlueckerVector, Nothing}`: The pluecker vector of the dual cell, required for linear and inverted linear dual cells.

# Throws
An error if the inputs are invalid.
"""
function check_dual_cell_inputs(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, cellType::Symbol)
    if cellType == :hypersurface
        @assert length(activeSupport) == 2 "active support must be a pair of indices for a hypersurface"
    elseif cellType == :linear || cellType == :inverted_linear
        # check that pluecker indices indexed by active support are loopless
        # assert that no row of the submatrix indexed by active support is zero
        @assert all([!iszero(v) for v in eachcol(ambientSupport[:, activeSupport]')]) "active support must be loopless"
    end
end

"""
    codim(m::DualCell)

Compute the codimension of a dual cell.

# Arguments
- `m::DualCell`: The dual cell.

# Returns
The codimension of the dual cell.
"""
function codim(m::DualCell)
    # compute the rank of the active support

    # take the submatrix of the ambient support indexed by active support
    indexedSupport = m.ambientSupport[m.activeSupport]

    # choose the first indexed element as a reference vector
    referenceVector = Matrix{Int}(vcat(m.ambientSupport[m.activeSupport[1], :])')
    # subtract this from every other row indexed by active support
    indexedSupport = indexedSupport .- referenceVector

    return rank(indexedSupport[2:end, :])
end

# handles printing of a dual cell
function Base.show(io::IO, m::DualCell)
    print(io, "Dual cell of type $(typeof(m).parameters[1]) supported on vertices $(m.activeSupport)")
end