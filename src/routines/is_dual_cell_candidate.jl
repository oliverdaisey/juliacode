include("../structs/dual_type.jl")
include("../structs/support.jl")

function is_dual_cell_candidate(S::Support{Hypersurface}, s::Vector{Int})
    return length(s) >= 2
end

function is_dual_cell_candidate(S::Support{<:Union{Linear, InvertedLinear}}, s::Vector{Int})
    return is_loopless(S, s)
end

function is_loopless(S::Support{<:Union{Linear, InvertedLinear}}, s::Vector{Int})
    Ss = S[s]
    return findfirst(iszero, [Ss[:,i] for i in 1:size(Ss,2)]) === nothing
end

function is_loopless(S::Support{<:Union{Linear, InvertedLinear}}, Ss::Vector{Vector{Int}})
    return findfirst(iszero, [Ss[:,i] for i in 1:size(Ss,2)]) === nothing
end

function is_dual_cell_candidate(S::Support{<:Union{Linear, InvertedLinear}}, Ss::Vector{Vector{Int}})
    return is_loopless(S, Ss)
end

function is_dual_cell_candidate(S::Support{Hypersurface}, Ss::Vector{Vector{Int}})
    return length(Ss) > 1
end

