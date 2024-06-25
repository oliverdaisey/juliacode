using Oscar

"""
    subdivision_of_points(points::Matrix{Int}, weights::Vector)

Given a matrix of points and a vector of weights, return the subdivision of the points. Temporary implmentation until the function with these arguments is correctly implemented in Oscar.
"""
function subdivision_of_points_workaround(points::Matrix{Int}, weights::Vector)
    indices = findall(!iszero, weights)

    relevantPoints = points[indices, :]
    relevantWeights = QQ.(weights[indices]; preserve_ordering=true)

    return subdivision_of_points(relevantPoints, relevantWeights)



end