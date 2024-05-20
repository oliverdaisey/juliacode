function random_lift(nu::TropicalSemiringMap{Kt,t,minOrMax}, a::Union{TropicalSemiringElem{minOrMax}, Nothing}=nothing) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    functionField = Oscar.valued_field(nu)
    coefficientField = base_ring(functionField)

    if isnothing(a)
        a = ZZ(rand(Int8))
    else
        a = ZZ(a; preserve_ordering=true)
    end
    randomLift = functionField(uniformizer(nu))^a

    if coefficientField==QQ
        return rand(-999:999)*randomLift
    else
        return rand(K)*randomLift
    end
end

"""
    random_lift(nu::TropicalSemiringMap{Kt, t, minOrMax}, f::AbstractAlgebra.Generic.MPoly{TropicalSemiringElem{minOrMax}}, parent::AbstractAlgebra.Generic.MPolyRing{AbstractAlgebra.Generic.RationalFunctionFieldElem{QQFieldElem, QQPolyRingElem}}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}

Randomly lift each coefficient of `f` and return the result in the correct polynomial ring.

# Inputs
- `nu::TropicalSemiringMap{Kt, t, minOrMax}`: A tropical semiring map.
- `f::AbstractAlgebra.Generic.MPoly{TropicalSemiringElem{minOrMax}}`: A tropical polynomial.
- `parent::AbstractAlgebra.Generic.MPolyRing{AbstractAlgebra.Generic.RationalFunctionFieldElem{QQFieldElem, QQPolyRingElem}}`: The ring where the lifted polynomials live.
"""
function random_lift(nu::TropicalSemiringMap{Kt, t, minOrMax}, f::AbstractAlgebra.Generic.MPoly{TropicalSemiringElem{minOrMax}}, parent::AbstractAlgebra.Generic.MPolyRing{AbstractAlgebra.Generic.RationalFunctionFieldElem{QQFieldElem, QQPolyRingElem}}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}

    # random lift each coefficient of f, return result in correct polynomial_ring
    return map_coefficients(c -> random_lift(nu, c), f, parent=parent)

end