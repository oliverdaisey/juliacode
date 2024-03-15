include("../structs/dual_type.jl")
include("../structs/support.jl")

function is_dual_cell_candidate(s::Vector{Int}, S::Support{Hypersurface})
    return length(s) >= 2
end

function is_dual_cell_candidate(s::Vector{Int}, S::Support{<:Union{Linear, InvertedLinear}})
    return is_loopless(s, S)
end

function is_loopless(s::Vector{Int}, S::Support{<:Union{Linear, InvertedLinear}})

    Ss = S[s]
    return findfirst(iszero, [Ss[:,i] for i in 1:size(Ss,2)]) === nothing
end

function is_loopless(Ss::Vector{Vector{Int}}, S::Support{<:Union{Linear, InvertedLinear}})
    return findfirst(iszero, [Ss[:,i] for i in 1:size(Ss,2)]) === nothing
end

function is_dual_cell_candidate(Ss::Vector{Vector{Int}}, S::Support{<:Union{Linear, InvertedLinear}})
    return is_loopless(Ss, S)
end

function is_dual_cell_candidate(Ss::Vector{Vector{Int}}, S::Support{Hypersurface})
    return length(Ss)>1
end