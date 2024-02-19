abstract type DualCellType end

struct DualCellHypersurface <: DualCellType end
struct DualCellLinear <: DualCellType end
struct DualCellInvertedLinear <: DualCellType end

struct DualCell{cellType<:DualCellType, minOrMax<:Union{typeof(min),typeof(max)}}
    ambientSupport :: Matrix{Int}
    activeSupport :: Vector{Int}

    function DualCell{cellType, minOrMax}(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}) where {cellType<:DualCellType, minOrMax<:Union{typeof(min),typeof(max)}}
        return new{cellType, minOrMax}(ambientSupport, activeSupport)
    end
end

r"""
    dual_cell(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, cellType::Symbol, ::typeof(min)=min, pluecker_vector::Union{TropicalPlueckerVector, Nothing}=nothing)

    Create a dual cell of the given type, using the min convention, with given ambient and active support.

    # Arguments
    - `ambientSupport::Matrix{Int}`: The ambient support of the dual cell.
    - `activeSupport::Vector{Int}`: The active support of the dual cell, corresponding to indices of columns of the ambient support.
    - `cellType::Symbol`: The type of the dual cell, must be one of :hypersurface, :linear, or :inverted_linear.
    - `::typeof(min)`: The min convention.
    - `pluecker_vector::Union{TropicalPlueckerVector, Nothing}`: The pluecker vector of the dual cell, required for linear and inverted linear dual cells.

    # Returns
    A dual cell of the given type, using the min convention, with given ambient and active support.
"""
function dual_cell(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, cellType::Symbol, ::typeof(min)=min, pluecker_vector::Union{TropicalPlueckerVector, Nothing}=nothing)
    check_dual_cell_inputs(ambientSupport, activeSupport, cellType, pluecker_vector)
    if cellType == :hypersurface
        return DualCell{DualCellHypersurface, typeof(min)}(ambientSupport, activeSupport)
    elseif cellType == :linear
        return DualCell{DualCellLinear, typeof(min)}(ambientSupport, activeSupport)
    elseif cellType == :inverted_linear
        return DualCell{DualCellInvertedLinear, typeof(min)}(ambientSupport, activeSupport)
    else
        throw(ArgumentError("cellType must be one of :hypersurface, :linear, or :inverted_linear"))
    end
end

r"""
    dual_cell(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, cellType::Symbol, ::typeof(max)=max, pluecker_vector::Union{TropicalPlueckerVector, Nothing}=nothing)

    Create a dual cell of the given type, using the max convention, with given ambient and active support.

    # Arguments
    - `ambientSupport::Matrix{Int}`: The ambient support of the dual cell.
    - `activeSupport::Vector{Int}`: The active support of the dual cell, corresponding to indices of columns of the ambient support.
    - `cellType::Symbol`: The type of the dual cell, must be one of :hypersurface, :linear, or :inverted_linear.
    - `::typeof(max)`: The max convention.
    - `pluecker_vector::Union{TropicalPlueckerVector, Nothing}`: The pluecker vector of the dual cell, required for linear and inverted linear dual cells.

    # Returns
    A dual cell of the given type, using the max convention, with given ambient and active support.
"""
function dual_cell(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, cellType::Symbol, ::typeof(max), pluecker_vector::Union{TropicalPlueckerVector, Nothing}=nothing)
    error("max tropical dual cells currently unsupported")
    check_dual_cell_inputs(ambientSupport, activeSupport, cellType, pluecker_vector)
    if cellType == :hypersurface
        return DualCell{DualCellHypersurface, typeof(max)}(ambientSupport, activeSupport)
    elseif cellType == :linear
        return DualCell{DualCellLinear, typeof(max)}(ambientSupport, activeSupport)
    elseif cellType == :inverted_linear
        return DualCell{DualCellInvertedLinear, typeof(max)}(ambientSupport, activeSupport)
    else
        throw(ArgumentError("cellType must be one of :hypersurface, :linear, or :inverted_linear"))
    end
end

r"""
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
function check_dual_cell_inputs(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, cellType::Symbol, pluecker_vector::Union{TropicalPlueckerVector, Nothing}=nothing)
    if cellType == :hypersurface
        @assert length(activeSupport) == 2 "active support must be a pair of indices for a hypersurface"
    elseif cellType == :linear || cellType == :inverted_linear
        if isnothing(pluecker_vector)
            error("pluecker vector required for linear or inverted linear dual cells")
        end
        # check that pluecker indices indexed by active support are loopless
        selected_indices = [pluecker_vector.pluecker_indices[i] for i in activeSupport]
        @assert length(unique(vcat((selected_indices)...))) == rank(pluecker_vector) "pluecker indices $(selected_indices) do not form a loopless set"
    end
end

r"""
    codim(m::DualCell)

    Compute the codimension of a dual cell.

    # Arguments
    - `m::DualCell`: The dual cell.

    # Returns
    The codimension of the dual cell.
"""
function codim(m::DualCell)
    # compute the rank of the active support
    return rank(m.ambientSupport[:,m.activeSupport])
end

# handles printing of a dual cell
function Base.show(io::IO, m::DualCell)
    print(io, "Dual cell of type $(typeof(m).parameters[1]) supported on vertices $(m.activeSupport)")
end