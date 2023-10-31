using Oscar
include("cayleyembedding.jl")

function polyhedral_lift(A::Matrix{Int}, ω::Vector{Int})
    # to do: add some checks on the input
    n = nrows(A) # number of variables
    m = ncols(A) # number of terms in equation
    
    point_list = []

    for i in 1:m
        pt = [A[j,i] for j in 1:n]
        push!(pt, ω[i])
        push!(point_list, pt)
    end
     
    upper_hull = convex_hull(point_list)

    # inequalities of the half space embedded in R^{n+1}

    equations_matrix = identity_matrix(QQ, n)
    vcat(equations_matrix, zeros(QQ, 1, n))
    hcat(equations_matrix, zeros(QQ, n+1, 1))

    return equations_matrix

end

A = [1 2; 3 4]
ω = [1,1]

polyhedral_lift(A, ω)