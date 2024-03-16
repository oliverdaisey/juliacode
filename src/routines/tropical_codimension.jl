include("../structs/dual_type.jl")
include("../structs/support.jl")

function tropical_codim(S::Support{Hypersurface})
    return 1
end

function tropical_codim(S::Support{<:Union{Linear, InvertedLinear}})
    return length(findall(x -> x != 0, S[1]))
end