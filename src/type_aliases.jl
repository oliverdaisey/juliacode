using Oscar

MixedCell = Vector{NTuple{2,Int}}
PointConfiguration = Vector{Matrix{QQFieldElem}}
GeneralizedMixedCell = Vector{Vector{Int}}
# alias for a tuple of tropical polynmomials
TropicalTuple{N} = NTuple{N, AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}}
TropicalPoly = AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}} # alias for a tropical polynomial
PluckerVector = Pair{Vector{Vector{Int}}, Vector{Int}}