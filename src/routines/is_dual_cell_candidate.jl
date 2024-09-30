
function is_dual_cell_candidate(S, s::Vector{Int})
    return length(s) >= 2
end

function is_dual_cell_candidate(S, Ss::Vector{Vector{Int}})
    return length(Ss) >= 2
end

# TODO: Write this with the new linear dual cell types
function is_dual_cell_candidate(S, s::Vector{Int})
    return is_loopless(S, s)
end

function is_dual_cell_candidate(S, Ss::Vector{Vector{Int}})
    return is_loopless(S, Ss)
end


function is_loopless(S, s::Vector{Int})
    return is_loopless(S, S[s])
end

function is_loopless(S, Ss::Vector{Vector{Int}})
    return findfirst(iszero, [[Ss[i][j] for i in 1:length(Ss)] for j in 1:length(Ss[1])]) === nothing
end
