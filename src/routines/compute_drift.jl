
# to do: use mixed cell cone of s to ensure we stay inside
function compute_drift(s::MixedCell, h::MixedPath, pointer_index::Int)
    initial_point = stable_intersection_point(s, h, pointer_index, QQFieldElem(0))
    final_point = stable_intersection_point(s, h, pointer_index, QQFieldElem(1//100))

    return final_point - initial_point
end