using Oscar
include("cayleyembedding.jl")

struct MixedCell
   indices::Vector{NTuple{2,Int}}
end



"""
Function to compute the bistellar flip as described in Anders' paper.

INPUTS:
- `A::Vector{Matrix{QQFieldElem}}`: a vector of point configurations in the same dimension
- `σ::MixedCell`: a mixed cell of the cayley embedding
- `c::Vector{QQFieldElem}`: a circuit defining a facet F of the mixed cell cone C_σ
- `A_mix::Vector{MixedCell}`: a vector of mixed cells sufficiently close to F
OUTPUTS:
- `B_mix::Vector{MixedCell}`: a vector of mixed cells sufficiently close to F after the bistellar flip

"""
function bistellar_flip(
    A::Vector{Matrix{QQFieldElem}}, σ::MixedCell, c::Vector{QQFieldElem}, A_mix::Vector{MixedCell}
) :: Vector{MixedCell} 
{

}
