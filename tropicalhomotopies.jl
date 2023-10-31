using Oscar
using("cayleyembedding.jl")

function polyhedral_lift(A::Matrix{Integer}, ω::Vector{QQ})
    n = nrows(A) # number of variables
    m = ncols(A) # number of terms in equation
    
    point_list = []

    for i in range(m)
        pt = [A[j,i] for j in range(n)]
        push!(pt, ω[i])
        push!(point_list, pt)
    end

    upper_hull = convex_hull(point_list)

    # inequalities of the half space embedded in R^{n+1}


end