
function perturb_dual_vector(S::DualSupport{Hypersurface}, c::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where {minOrMax <: Union{typeof(min), typeof(max)}}
    epsilon = QQ(1//1000000)
    perturbation = epsilon.*rand(-99:99, length(c))

    return c .* parent(first(c)).(perturbation)
end

function perturb_dual_vector(S::DualSupport{<:Union{Linear, InvertedLinear}}, c::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where {minOrMax <: Union{typeof(min), typeof(max)}}

    epsilon = QQ(1//1000000)
    n = tropical_ambient_dim(S)
    tropicalPerturbation = epsilon.*rand(-99:99, n)

    @assert length(c) == length(S) "The length of the dual vector must be equal to the length of the ambient support."

    perturbation = zeros(QQ, length(c))
    for (i, c) in enumerate(c) 
        Si = S[i]
        perturbation[i] = sum(tropicalPerturbation[j] for j in 1:n if Si[j] != 0)
    end

    return c .* parent(first(c)).(perturbation)
end

