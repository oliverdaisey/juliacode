
mutable struct LazyLinearDualCell{cellType<:Union{Linear, InvertedLinear},minOrMax<:Union{typeof(min),typeof(max)}}

    ambientDualSupport::LazyTropicalPlueckerVector{minOrMax}
    activeIndices::Vector{Vector{Int}}

    function LazyLinearDualCell{cellType,minOrMax}(ambientDualSupport::LazyTropicalPlueckerVector{minOrMax}, activeIndices::Vector{Vector{Int}}) where {cellType<:Union{Linear, InvertedLinear},minOrMax<:Union{typeof(min),typeof(max)}}
        return new{cellType,minOrMax}(ambientDualSupport, activeIndices)
    end
end