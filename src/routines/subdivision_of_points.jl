using Oscar

function subdivision_of_points_workaround(points::Matrix{Int}, weights::Vector)
    indices = findall(!iszero, weights)

    relevantPoints = points[indices, :]
    relevantWeights = QQ.(weights[indices]; preserve_ordering=true)

    return subdivision_of_points(relevantPoints, relevantWeights)



end