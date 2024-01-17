using Oscar

"""
ω is generic means it is not in the codimension 2 skeleton of the secondary fan of
cayley_embedding
"""

# type aliases
const MixedCell = Vector{NTuple{2,Int}}
const PointConfiguration = Vector{Matrix{QQFieldElem}}

struct CayleyCircuit
    # coefficients defining the corresponding facet
    coefficients::Vector{QQFieldElem}

    # indices of the columns supporting this circuit
    support::Vector{Int}

    # configuqation and column index in cayley embedding giving rise to this circuit
    i::Int
    gamma::Int
end

"""
Utility to remove an item from a vector
"""
function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

"""
Function to compute the bistellar flip as described in Anders' paper.

INPUTS:
- `A::Vector{PointConfiguration}`: a vector of point configurations in the same dimension
- `σ::MixedCell`: a mixed cell of the cayley embedding. We encode it as a choice of two indices from each
column of the cayley embedding in each dimension.
- `c::Vector{QQFieldElem}`: a circuit defining a facet F of the mixed cell cone C_σ
- `A_mix::Vector{MixedCell}`: a vector of mixed cells sufficiently close to F
OUTPUTS:
- `B_mix::Vector{MixedCell}`: a vector of mixed cells sufficiently close to F after the bistellar flip

"""
function bistellar_flip(
    A::Vector{PointConfiguration}, σ::MixedCell, c::CayleyCircuit, A_mix::Vector{MixedCell}
)::Vector{MixedCell}

    # step 1
    B_mix = Vector{MixedCell}()

    # step 2 is handled by the caller
    i = c.I
    gamma = c.gamma

    # step 3
    α, β = σ[i]

    # step 4
    for τ in A_mix
        if gives_rise(τ, c)
            mixed_cell_cayley_indices = reduce(vcat, σ)
            if c.coefficients[α] > 0
                updated_indices = deepcopy(mixed_cell_cayley_indices)
                push!(updated_indices, gamma)
                remove!(updated_indices, α)
                if updated_indices ∉ B_mix
                    push!(B_mix, updated_indices)
                end 
            end
            # note we can update B_mix twice in one step!
            if c.coefficients[β] > 0
                updated_indices = deepcopy(mixed_cell_cayley_indices)
                push!(updated_indices, gamma)
                remove!(updated_indices, β)
                if updated_indices ∉ B_mix
                    push!(B_mix, updated_indices)
                end 
            end
        end
    end

    # step 5
    return B_mix
    
end

"""
Decides whether a mixed cell τ gives rise to a circuit
"""
function gives_rise(τ::MixedCell, c::CayleyCircuit)::Bool
    
end

