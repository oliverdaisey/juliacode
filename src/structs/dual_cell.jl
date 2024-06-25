
mutable struct DualCell{cellType<:DualType, minOrMax<:Union{typeof(min), typeof(max)}}

    dualWeight::DualWeight
    activeSupport::Vector{Vector{Int}}

end


"""
    dual_cell(ambientDualSupport::Matrix{Int}, activeIndices::Vector{Int}, dualWeight::Vector{Oscar.TropicalSemiringElem{minOrMax}}, cellType::Symbol, ::typeof(min)=min)

Create a dual cell of the given type, using the min convention, with given ambient and active support. This is the entry point for constructing all dual cells.

# Arguments
- `ambientDualSupport::Matrix{Int}`: The ambient support of the dual cell, with columns as points.
- `activeIndices::Vector{Int}`: The active support of the dual cell, corresponding to indices of columns of the ambient support.
- `dualWeight::Vector{Oscar.TropicalSemiringElem{minOrMax}}`: The dual vector giving rise to this dual cell.
- `cellType::Symbol`: The type of the dual cell, must be one of :hypersurface, :linear, or :inverted_linear.
- `::typeof(min)`: The min convention.

# Returns
A dual cell of the given type, using the min convention, with given ambient and active support.
"""
function dual_cell(activeSupport::Vector{Vector{Int}}, dualWeight::DualWeight{Hypersurface, typeof(min)}, ::typeof(min)=min)
    check_dual_cell_inputs(Hypersurface, activeSupport)
    return DualCell{Hypersurface, typeof(min)}(dualWeight, activeSupport)
end

function dual_cell(activeSupport::Vector{Vector{Int}}, dualWeight::DualWeight{LinearType, typeof(min)}, ::typeof(min)=min) where (LinearType<:Union{Linear, InvertedLinear})
    check_dual_cell_inputs(LinearType, activeSupport)
    if LinearType == Linear
        return DualCell{Linear, typeof(min)}(dualWeight, activeSupport)
    else
        return DualCell{InvertedLinear, typeof(min)}(dualWeight, activeSupport)
    end
end

"""
    dual_cell(dualWeight::DualWeight{Hypersurface, typeof(min)}, ::typeof(min)=min)

Assume that all points are active.
"""
function dual_cell(dualWeight::DualWeight{Hypersurface, typeof(min)}, ::typeof(min)=min)
    # construct active support, then call the other constructor
    polynomial = data(dualWeight)
    return dual_cell(generate_support(polynomial), dualWeight)

end

"""
    dual_cell(polynomial::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}, ::typeof(min)=min)

Create a hypersurface dual cell from a tropical polynomial. All points are active.
"""
function dual_cell(polynomial::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}, ::typeof(min)=min)
    return dual_cell(generate_support(polynomial), dual_weight(polynomial))
end



# linear dual cell constructor given realisation matrix, tropical semiring map nu, and a vector of tropical semiring elems inside the tropical linear space
function dual_cell(dualType::LinearType, realisation::MatElem{Kelem}, nu::TropicalSemiringMap{K,t,minOrMax}, interiorPoint::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where {Kelem, K,t,minOrMax, LinearType<:Union{typeof(Linear), typeof(InvertedLinear)}}

    @assert Kelem == elem_type(K) "The realisation must have entries of the correct type"
    dualWeight = dual_weight(dualType, realisation, nu)

    # work out active support based on the given point
    flats = [findall(x -> x == val, interiorPoint) for val in unique(interiorPoint)]

    return DualCell{dualType, typeof(min)}(dualWeight, activeSupport)
end



function check_dual_cell_inputs(cellType, activeSupport::Vector{Vector{Int}}) 

    # We just return for now, until I write this function properly.
    return

    if cellType == :hypersurface
        @assert length(activeIndices) == 2 "active support must be a pair of indices for a hypersurface"
    elseif cellType == :linear || cellType == :inverted_linear
        # check that pluecker indices indexed by active support are loopless
        # assert that no row of the submatrix indexed by active support is zero
        @assert all([!iszero(v) for v in eachcol(ambientDualSupport[:, activeIndices]')]) "active support must be loopless"
    end
end

"""
    tropical_codim(m::DualCell)

Computes the codimension of the tropical polyhedron dual to `m`.
"""
function tropical_codim(m::DualCell)
    return Oscar.dim(convex_hull(active_support(m)))
end 

import Oscar.type

function type(m::DualCell)
    return typeof(m).parameters[1]
end

# handles printing of a dual cell
function Base.show(io::IO, m::DualCell)
    print(io, "$(type(m)) dual cell supported on $(length(active_support(m))) points")
end

function ambient_dim(m::DualCell)
    return length(first(active_support(m)))
end

function rank(m::DualCell{LinearType, typeof(min)}) where LinearType<:Union{typeof(Linear), typeof(InvertedLinear)}
    # return number of non zero elements of first(active_support(m))
    return sum([!iszero(v) for v in first(active_support(m))])
end

function active_support(m::DualCell)
    return m.activeSupport
end

function dual_weight(m::DualCell)
    return m.dualWeight
end

function convention(m::DualCell{<:DualType, typeof(min)})
    return min
end

function convention(m::DualCell{<:DualType, typeof(max)})
    return max
end

function center(m::DualCell)
    return QQ.([sum(active_support(m), dims=1)...]) / length(active_support(m))
end

import Oscar.convex_hull
function convex_hull(m::DualCell)
    return convex_hull(active_support(m))
end

function dual_facets(m::DualCell)
    cellHull = convex_hull(active_support(m))
    properDualFaces = DualCell[]
    for face in faces(cellHull, Oscar.dim(cellHull) - 1)
        faceVertices = [Int.(v) for v in vertices(face)]
        # println(faceVertices)
        # println("is dual cell candidate: ", is_dual_cell_candidate(ambient_support(m), faceVertices))
        if is_dual_cell_candidate(ambient_support(m), faceVertices)
            push!(properDualFaces, dual_cell(ambient_support(m), faceVertices, dual_weight(m)))
        end
    end
        
    return properDualFaces
end

"""
    dual_cells(S::DualSupport{<:DualType}, c::Vector{Oscar.TropicalSemiringElem})

Create a vector of all the dual cells from a support (of exponent vectors) and a lift.

"""
function dual_cells(S, c::Vector{<:Oscar.TropicalSemiringElem})

    dualSubdivision = subdivision_of_points_workaround(points(S), c)
    polyhedralComplex = polyhedral_complex(Oscar.pm_object(dualSubdivision).POLYHEDRAL_COMPLEX)

    dualCells = []
    for i in tropical_codim(S):dim(polyhedralComplex)
        for p in polyhedra_of_dim(polyhedralComplex, i)
            s = [Int.(v) for v in vertices(p)]
            if is_dual_cell_candidate(S, s)
                push!(dualCells, dual_cell(S, s, c))
            end
        end
    end

    return dualCells
end

function generate_support(f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}})
    return collect(exponents(f))
end