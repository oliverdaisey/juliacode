using Oscar

"""
    generate_support(f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}})

Generate the support of a tropical polynomial.

# Arguments
- `f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}`: The tropical polynomial.

# Returns
- The support of the tropical polynomial.
"""
function generate_support(f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}})::Matrix{Int64}
    return Matrix{Int}((f.exps)')
end