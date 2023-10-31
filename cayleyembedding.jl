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
    @assert length(unique(length.(vcat(P1,P2)))) == 1 "points in P1 and P2 must have same dimension"

    P1padded = [ vcat(p1, [1,0]) for p1 in P1 ]
    P2padded = [ vcat(p2, [0,1]) for p2 in P2 ]
    return vcat(P1padded, P2padded)

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

P1 = [[0,0], [1,0], [0,1], [1,1]]
P2 = [[0,0], [1,0], [2,1]]
w1 = [0,-1,-1,-2]
w2 = [0,0,0]
minkowski_subdivision, minkowski_labels = minkowski_sum_regular_subdivision(P1, P2, w1, w2)
println(is_minkowski_label_mixed.(values(minkowski_labels)))

P1 = [[0,0],[1,0],[0,1]]
P2 = [[0,0],[1,0],[0,1]]
w1 = [0,0,0]
w2 = [0,-1,0]
minkowski_subdivision, minkowski_labels = minkowski_sum_regular_subdivision(P1, P2, w1, w2)
println(is_minkowski_label_mixed.(values(minkowski_labels)))
