"""
    ray_intersects_cone(C::Cone, u::Vector{QQFieldElem}, v::Vector{QQFieldElem})

Determines the value of ``t`` at which the ray ``u + tv`` with ``t in [0, 1]`` first intersects the boundary of the cone ``C``.

INPUTS:
- ``C::Cone``: cone
- ``u::Vector{QQFieldElem}``: start point of the ray, required to be in C
- ``v::Vector{QQFieldElem}``: direction of the ray, required to be nonzero

OUTPUTS:
- ``t::QQFieldElem``: parameter value at which the ray first intersects the boundary of the cone, or ``nothing`` if the ray does not intersect the cone.

# Examples
```jldoctest
julia> C = cone_from_inequalities([-1 0 0; 0 -1 0; 0 0 -1])
julia> u = Vector{QQFieldElem}([1, 1, 1])
julia> v = Vector{QQFieldElem}([-1, 0, 0])
julia> ray_intersects_cone(C, u, v)
1
julia> v[1] = 1
julia> ray_intersects_cone(C, u, v)
nothing
```
"""
function ray_intersects_cone(C::Polyhedron, u::Vector{QQFieldElem}, v::Vector{QQFieldElem})

    @req u in C "u must be in the cone C"
    @req !iszero(v) "v must be nonzero"
    
    intersectionTimes = []
    for F in facets(C)
        v_dot = dot(F.a , v)
        if v_dot == 0
            # this means they are parallel, so we can't intersect
            continue
        end
        t = dot(-F.a , u) / v_dot
        if t > 0
            push!(intersectionTimes, t)
        end
    end
    
    if isempty(intersectionTimes)
        # nothing satisfied the criteria, so ray does not intersect
        return nothing
    else
        # intersection time is precisely the minimum of the intersection points
        return minimum(intersectionTimes)
    end


end