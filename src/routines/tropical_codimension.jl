
function tropical_codim(S::DualSupport{Hypersurface})
    return 1
end

function tropical_codim(S::DualSupport{<:Union{Linear, InvertedLinear}})
    return length(findall(x -> x != 0, S[1]))
end