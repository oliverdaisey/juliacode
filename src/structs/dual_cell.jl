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

function dual_cell(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, cellType::Symbol, ::typeof(min)=min)
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

function dual_cell(ambientSupport::Matrix{Int}, activeSupport::Vector{Int}, cellType::Symbol, ::typeof(max))
    error("max tropical dual cells currently unsupported")
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

# todo: dual cell printing
