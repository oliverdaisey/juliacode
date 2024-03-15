include("../structs/dual_type.jl")
include("../structs/support.jl")

function is_dual_cell_candidate(s::Vector{Int}, S::Support{<:DualType})
    type = typeof(S).parameters[1]
    if type == Hypersurface
        return length(s) >= 2
    elseif type == Linear || S.parameters == InvertedLinear
        return is_loopless(s, S)
    else
        throw(ArgumentError("Support type must be one of :Hypersurface, :Linear, or :InvertedLinear"))
    end
end

function is_loopless(s::Vector{Int}, S::Support{<:Union{Linear, InvertedLinear}})
    throw(ArgumentError("is_loopless not yet implemented"))
end