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
