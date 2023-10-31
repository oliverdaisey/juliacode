using Oscar
include("cayleyembedding.jl")

"""
Given a matrix of exponent vectors and coefficient vector ω, returns the polyhedral lift of A induced by ω.
# Arguments:
- `A::Matrix{Int}`: a matrix of exponent vectors
- `ω::Vector{Int}`: a vector of coefficients
# Returns:
a polyhedron in R^{n+1} that is the polyhedral lift of A induced by ω.
"""
function polyhedral_lift(A::Matrix{Int}, ω::Vector{Int})
    
    # check inputs
    @req ncols(A) == length(ω) "number of columns in A must match length of ω"

    n = nrows(A) # number of variables
    m = ncols(A) # number of terms in equation
    
    point_list = Vector{Vector{Int}}()

    for i in 1:m
        pt = [A[j,i] for j in 1:n]
        push!(pt, ω[i])
        push!(point_list, pt)
    end
    
    upper_hull = convex_hull(point_list)

    # equations and inequalities of the half space embedded in R^{n+1}

    equations_matrix = identity_matrix(QQ, n)
    equations_matrix = vcat(equations_matrix, matrix(QQ, zeros(QQ, 1, n)))
    equations_matrix = hcat(equations_matrix, matrix(QQ, zeros(QQ, n+1, 1)))

    inequalities_matrix = zero_matrix(QQ, n+1, n+1)
    inequalities_matrix[n+1, n+1] = 1

    println(inequalities_matrix)

    # b = zero_matrix(QQ, n+1, 1)
    # b = [b...] # convert b to a form that the polymake interface likes
    b = zeros(QQ, n+1)

    half_line = polyhedron((inequalities_matrix, b), (equations_matrix, b))

    return upper_hull + half_line

end