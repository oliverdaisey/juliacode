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

"""
Minkowski sum of a vector of point configurations, with optional weights.
"""
function minkowski_sum_points(M::Vector{Matrix{QQFieldElem}}, W::T=nothing) where {T<:Union{Nothing, Vector{Vector{Int}}}}

    @req length(unique(ncols.(M))) == 1 "points in M must have same dimension"
    @req isnothing(W) || length(W) == length(M) "number of weight vectors must match number of point configurations"
    @req isnothing(W) || unique([length(W[i]) == nrows(M[i]) for i in 1:length(W)]) == [true] "number of weights in weight vector must match number of points in point configuration"

    ms = Vector{QQFieldElem}[]
    m = length(M) # number of point configurations

    for I in Iterators.product([1:nrows(M[i]) for i in 1:m]...)
        vec = [M[i][I[i],:] for i in 1:m]
        push!(ms, sum(vec))
    end

    uniquems = unique(ms)

    if T != Nothing
        ws = Int[]
        for I in Iterators.product([1:length(W[i]) for i in 1:length(W)]...)
            push!(ws, sum([W[i][I[i]] for i in 1:length(W)]))
        end
        uniquews = [reduce(min, ws[findall(x -> x == pt, ms)]) for pt in uniquems]
        return uniquems, uniquews
    end

    return uniquems

end


function minkowski_sum_points(P1::Vector{Vector{Int}}, P2::Vector{Vector{Int}}, w1::T=nothing, w2::T=nothing) where {T<:Union{Nothing, Vector{Int}}}

    @assert length(unique(length.(vcat(P1,P2)))) == 1 "points in P1 and P2 must have same dimension"
    @assert isnothing(w1) || length(w1) == length(P1) "weights must match number of points in configuration"
    @assert isnothing(w2) || length(w2) == length(P2) "weights must match number of points in configuration"

    ms = Vector{Int}[]
    for v in P1
        for w in P2
            push!(ms, v+w)
        end
    end
    uniquems = unique(ms)

    if T != Nothing
        ws = Int[]
        for w1i in w1
            for w2j in w2
                push!(ws, w1i + w2j)
            end
        end
        uniquews = [reduce(min, ws[findall(x -> x == pt, ms)]) for pt in uniquems]
        return uniquems, uniquews
    end

    return uniquems

end

function minkowski_sum_regular_subdivision(M::Vector{Matrix{QQFieldElem}}, W::Vector{Vector{Int}})
    
    cayley_points = cayley_embedding(M)
    cayley_subdivision = subdivision_of_points(matrix(QQ, cayley_points), vcat(W...))

    minkowski_points = minkowski_sum_points(M)
    minkowski_cells = Vector{Int}[]
    minkowski_labels = Vector{Int}[]
    n = ncols(M[1]) # dimension of points
    m = length(M) # number of point configurations

    for sigma in maximal_cells(cayley_subdivision)
        sigma_points = cayley_points[sigma]
        println(sigma_points)
        sigma_points_flattened = []
        for i in 1:m
            if sigma_points[n+i]==1
                push!(sigma_points_flattened, sigma_points[1:n])
            end
        end
        minkowski_cell_points = Vector{Int}[]
        for I in Iterators.product([1:length(sigma_points_flattened[i]) for i in 1:m]...)
            vec = [sigma_points_flattened[i][I[i]] for i in 1:m]
            push!(minkowski_cell_points, sum(vec))
        end

        minkowski_cell = Int[]
        for minkowski_cell_point in minkowski_cell_points
            push!(minkowski_cell, findfirst(isequal(minkowski_cell_point), minkowski_points))
        end

        minkowski_label = [dim(convex_hull(sigma_points_flattened[i])) for i in 1:m]

        push!(minkowski_cells, minkowski_cell)
        push!(minkowski_labels, minkowski_label)
    end

    minkowski_subdivision = subdivision_of_points(matrix(QQ,minkowski_points), IncidenceMatrix(minkowski_cells))
    minkowski_labels = Dict(zip(minkowski_cells,minkowski_labels))

    return minkowski_subdivision, minkowski_labels


end

function minkowski_sum_regular_subdivision(P1::Vector{Vector{Int}}, P2::Vector{Vector{Int}}, w1::Vector{Int}, w2::Vector{Int})

    ###
    # Construct subdivision of the Cayley embedding
    ###
    cayley_points = cayley_embedding(P1, P2)
    cayley_subdivision = subdivision_of_points(matrix(QQ,cayley_points), vcat(w1, w2))


    minkowski_points = minkowski_sum_points(P1,P2)
    minkowski_cells = Vector{Int}[]
    minkowski_labels = Vector{Int}[]
    n = length(first(P1))

    for sigma in maximal_cells(cayley_subdivision)
        sigma_points = cayley_points[sigma]
        sigma_points_flattened = [ [sigma_point[1:n] for sigma_point in sigma_points if sigma_point[n+i]==1 ] for i in 1:2]
        minkowski_cell_points = Vector{Int}[]
        for v in sigma_points_flattened[1]
            for w in sigma_points_flattened[2]
                vw = v+w
                if vw ∉ minkowski_cell_points
                    push!(minkowski_cell_points, v+w)
                end
            end
        end

        minkowski_cell = Int[]
        for vw in minkowski_cell_points
            push!(minkowski_cell, findfirst(isequal(vw), minkowski_points))
        end

        minkowski_label = [dim(convex_hull(sigma_points_flattened[i])) for i in 1:2]

        push!(minkowski_cells, minkowski_cell)
        push!(minkowski_labels, minkowski_label)
    end

    minkowski_subdivision = subdivision_of_points(matrix(QQ,minkowski_points), IncidenceMatrix(minkowski_cells))
    minkowski_labels = Dict(zip(minkowski_cells,minkowski_labels))

    return minkowski_subdivision, minkowski_labels
end


function is_minkowski_label_mixed(minkowski_label::Vector{Int})
    return findfirst(iszero,minkowski_label)===nothing
end
