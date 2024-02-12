using Base: BitSignedSmall, length_continued
using Oscar

# cayley embedding code
# given two point configurations P1 and P2, find their cayley embedding

"""
This function computes the Cayley embedding of two point configurations.
# Arguments:
- `p1::Vector{Vector{Int}}`: a point configuration
- `p2::Vector{Vector{Int}}`: a point configuration

# Returns:
the point configuration of the Cayley embedding.
"""
function cayley_embedding(P1::Vector{Vector{Int}}, P2::Vector{Vector{Int}})
    @req length(unique(length.(vcat(P1,P2)))) == 1 "points in P1 and P2 must have same dimension"

    P1padded = [ vcat(p1, [1,0]) for p1 in P1 ]
    P2padded = [ vcat(p2, [0,1]) for p2 in P2 ]
    return vcat(P1padded, P2padded)

end

"""
Cayley embedding given a vector of point configurations.
"""
function cayley_embedding(N::Vector{QQMatrix})

    M = copy(N)

    @req length(unique(nrows.(M))) == 1 "points in M must have same dimension"

    m = length(M) # number of point configurations
    

    for i in 1:m
        k = ncols(M[i]) # number of points in configuration i
        rows_to_append = matrix(QQ, zeros(QQ, m, k))
        rows_to_append[i,:] = ones(Int, k)
        M[i] = vcat(M[i], rows_to_append)
    end

    return reduce(hcat, M)

end