using Oscar

"""
Cayley embedding given a vector of point configurations. The points should be given as row vectors.
"""
function cayley_embedding(N::Vector{QQMatrix})

    M = copy(N)

    @req length(unique(ncols.(M))) == 1 "points in M must have same dimension"

    m = length(M) # number of point configurations
    

    for i in 1:m
        k = nrows(M[i]) # number of points in configuration i
        cols_to_append = matrix(QQ, zeros(QQ, k, m))
        cols_to_append[:,i] = ones(Int, k)
        M[i] = hcat(M[i], cols_to_append)
    end

    return reduce(vcat, M)

end