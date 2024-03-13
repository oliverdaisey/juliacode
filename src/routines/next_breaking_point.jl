include("../structs/mixed_cell_cone.jl")
include("ray_intersects_cone.jl")

"""
    next_breaking_point(s::MixedCell, h::MixedPath, pointer_index::Int, fraction::QQFieldElem)

    Given a mixed cell `s`, a mixed path `h`, a pointer index `pointer_index`, and a fraction `fraction` along the path,
    this function returns the next breaking point of the path relative to the given pointer index and fraction.
"""
function next_breaking_point(s::MixedCell, h::MixedPath, pointer_index::Int, fraction::QQFieldElem)
    
    dual_path_pointers = h.pointers[pointer_index]
    n = length(dual_path_pointers)

    # get the lift at this time
    lift = vcat([lift_from_node_and_fraction(h.dualPaths[i], dual_path_pointers[i], fraction) for i in 1:n]...)

    # get direction path is travelling in dual space
    next_dual_path_pointers = h.pointers[pointer_index + 1]
    direction = vcat([h.dualPaths[i].nodes[next_dual_path_pointers[i]] for i in 1:n]...) ./ vcat([h.dualPaths[i].nodes[dual_path_pointers[i]] for i in 1:n]...)

    C_s = mixed_cell_cone(s)

    # starting at `lift` and moving in the direction `direction`, when do we hit facet of C_s?
    # we want to find the smallest t such that lift + t*direction is on the facet of C_s
    t = ray_intersects_cone(mixed_cell_cone_to_polyhedron(C_s), lift, [x.data for x in direction])

    if isnothing(t)
        return nothing
    end

    return t - fraction

end