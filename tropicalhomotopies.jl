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

function normal_complex(A::Matrix{Int}, ω::Vector{Int})

    n = nrows(A)
    fan = normal_fan(polyhedral_lift(A, ω))
    
    # compute the RHS of the wedge
    equations_matrix = zeros(QQ, n+1)
    equations_matrix[n+1] = 1

    plane = polyhedron((zeros(QQ, n+1), [0]), (equations_matrix, 1))
    
    plane_to_complex = polyhedral_complex(
        IncidenceMatrix([[1]]),
        equations_matrix,
        Vector{Int}(),
        lineality_space(plane)
    )

    return common_refinement(polyhedral_complex_from_fan(fan), plane_to_complex)
    # to do: return the projection forgetting the last coordinate

end

function polyhedral_complex_from_fan(Σ::PolyhedralFan)

    rayvectors = [convert(Vector{QQFieldElem}, ray) for ray in rays(Σ)]
    
    return polyhedral_complex(
        maximal_cones(IncidenceMatrix, Σ),
        rayvectors, 
        [i for i in 1:length(rayvectors)],
        lineality_space(Σ))
    
end

A = [0 0 1 1; 0 2 0 1]
ω = [0, 0, 0, -2]

refinement = normal_complex(A, ω)