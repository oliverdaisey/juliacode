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

    return polyhedral_complex_from_fan(fan)

    return common_refinement(polyhedral_complex_from_fan(fan), plane_to_complex)
    # to do: return the projection forgetting the last coordinate
    # WARNING: functions involved here do not work as expected!!

end

function polyhedral_complex_from_fan(Σ::PolyhedralFan)

    rayvectors = [convert(Vector{QQFieldElem}, ray) for ray in rays(Σ)]
    
    return polyhedral_complex(
        maximal_cones(IncidenceMatrix, Σ),
        rayvectors, 
        [i for i in 1:length(rayvectors)],
        lineality_space(Σ))
    
end

function polyhedral_complex_from_polyhedra(Sigma::Vector{<:Polyhedron})
    ###
    # Construct matrix and incidence matrix of vertices and rays
    ###
    SigmaVertexIncidences = Vector{Int}[]
    SigmaVertices = Vector{QQFieldElem}[]
    SigmaRayIncidences = Vector{Int}[]
    SigmaRays = Vector{QQFieldElem}[]
    for sigma in Sigma
        sigmaVertexIncidence = Int[]
        for vertex in vertices(sigma)
            i = findfirst(isequal(vertex),SigmaVertices)
            if i === nothing
                push!(SigmaVertices,vertex)
                push!(sigmaVertexIncidence,length(SigmaVertices))
            else
                push!(sigmaVertexIncidence,i)
            end
        end
        push!(SigmaVertexIncidences,sigmaVertexIncidence)

        sigmaRayIncidence = Int[]
        for ray in rays(sigma)
            i = findfirst(isequal(ray),SigmaRays)
            if i === nothing
                push!(SigmaRays,ray)
                push!(sigmaRayIncidence,length(SigmaRays))
            else
                push!(sigmaRayIncidence,i)
            end
        end
        push!(SigmaRayIncidences,sigmaRayIncidence)
    end

    ###
    # Concatenate vertically matrixes of vertices and rays,
    # shift incidence matrix of rays and concatenate it horizontally to incicende matrix of vertices,
    # dehomogenize generators of lineality space
    ###
    SigmaVerticesAndRays = matrix(QQ,vcat(SigmaVertices,SigmaRays))
    SigmaRayIncidences = (x -> x .+length(SigmaVertices)).(SigmaRayIncidences)
    SigmaVertexAndRayIncidences = IncidenceMatrix([vcat(iv,ir) for (iv,ir) in zip(SigmaVertexIncidences,SigmaRayIncidences)])

    ###
    # Dehomogenize lineality space
    ###
    SigmaLineality = matrix(QQ,lineality_space(first(Sigma)))

    return polyhedral_complex(SigmaVertexAndRayIncidences,
                              SigmaVerticesAndRays,
                              collect(length(SigmaVertices)+1:nrows(SigmaVerticesAndRays)),
                              SigmaLineality)
end

# driver code
A = [0 0 1 1; 0 2 0 1]
ω = [0, 0, 0, -2]

refinement = normal_complex(A, ω)
maximal_polyhedra(refinement)