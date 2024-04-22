function random_lift(nu::TropicalSemiringMap{Kt,t,minOrMax}, a::TropicalSemiringElem{minOrMax}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    functionField = Oscar.valued_field(nu)
    coefficientField = base_ring(functionField)

    a = QQ(a; preserve_ordering=true)
    @assert isone(denominator(a)) "only tropical integers support"
    randomLift = functionField(uniformizer(nu))^ZZ(a)

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
    return map_coefficients(x -> random_lift(nu, x), f, parent=parent)

end