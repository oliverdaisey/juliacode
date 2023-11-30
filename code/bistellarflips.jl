using Oscar

# type aliases
const MixedCell = Vector{NTuple{2,Int}}
const PointConfiguration = Vector{Matrix{QQFieldElem}}


"""
Function to compute the bistellar flip as described in Anders' paper.

INPUTS:
- `A::Vector{PointConfiguration}`: a vector of point configurations in the same dimension
- `σ::MixedCell`: a mixed cell of the cayley embedding
- `c::Vector{QQFieldElem}`: a circuit defining a facet F of the mixed cell cone C_σ
- `A_mix::Vector{MixedCell}`: a vector of mixed cells sufficiently close to F
OUTPUTS:
- `B_mix::Vector{MixedCell}`: a vector of mixed cells sufficiently close to F after the bistellar flip

"""
function bistellar_flip(
    A::Vector{PointConfiguration}, σ::MixedCell, c::Vector{QQFieldElem}, A_mix::Vector{MixedCell}
)::Vector{MixedCell}
    B_mix = Vector{MixedCell}()

    return B_mix
    
end

