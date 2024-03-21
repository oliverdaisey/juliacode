"""
    DualCell{cellType<:DualType,minOrMax<:Union{typeof(min),typeof(max)}}

A dual cell of the given type, using the given convention.
"""
mutable struct DualCell{cellType<:DualType,minOrMax<:Union{typeof(min),typeof(max)}}

    ambientDualSupport::DualSupport{cellType}
    activeIndices::Vector{Int}
    dualWeight::Vector{Oscar.TropicalSemiringElem{minOrMax}}

    function DualCell{cellType,minOrMax}(ambientDualSupport::DualSupport{cellType}, activeIndices::Vector{Int}, dualWeight::Vector{Oscar.TropicalSemiringElem{minOrMax}}) where {cellType<:DualType,minOrMax<:Union{typeof(min),typeof(max)}}
        return new{cellType,minOrMax}(ambientDualSupport, activeIndices, dualWeight)
    end
end

"""
    dual_cell(ambientDualSupport::Matrix{Int}, activeIndices::Vector{Int}, dualWeight::Vector{Oscar.TropicalSemiringElem{minOrMax}}, cellType::Symbol, ::typeof(min)=min)

Create a dual cell of the given type, using the min convention, with given ambient and active support.

# Arguments
- `ambientDualSupport::Matrix{Int}`: The ambient support of the dual cell, with columns as points.
- `activeIndices::Vector{Int}`: The active support of the dual cell, corresponding to indices of columns of the ambient support.
- `dualWeight::Vector{Oscar.TropicalSemiringElem{minOrMax}}`: The dual vector giving rise to this dual cell.
- `cellType::Symbol`: The type of the dual cell, must be one of :hypersurface, :linear, or :inverted_linear.
- `::typeof(min)`: The min convention.

# Returns
A dual cell of the given type, using the min convention, with given ambient and active support.
"""
function dual_cell(ambientDualSupport::Matrix{Int}, activeIndices::Vector{Int}, dualWeight::Vector{Oscar.TropicalSemiringElem{typeof(min)}}, cellType::Symbol, ::typeof(min)=min)
    check_dual_cell_inputs(ambientDualSupport, activeIndices, cellType)
    if cellType == :hypersurface
        return DualCell{Hypersurface,typeof(min)}(ambientDualSupport, activeIndices)
    elseif cellType == :linear
        return DualCell{Linear,typeof(min)}(ambientDualSupport, activeIndices)
    elseif cellType == :inverted_linear
        return DualCell{InvertedLinear,typeof(min)}(ambientDualSupport, activeIndices)
    else
        throw(ArgumentError("cellType must be one of :hypersurface, :linear, or :inverted_linear"))
    end
end

function dual_cell(S::DualSupport{cellType}, s::Vector{Vector{Int}}, dualWeight::Vector{<:Oscar.TropicalSemiringElem}) where cellType<:DualType
    activeIndices = Int[]
    for alpha in s
        activeIndex = findfirst(i -> S[i] == alpha, 1:size(S, 1))
        @assert !isnothing(activeIndex) "Active support must be a subset of ambient support"
        push!(activeIndices, activeIndex)
    end

    return DualCell{cellType, typeof(min)}(S, activeIndices, dualWeight)
end


"""
    dual_cell(ambientDualSupport::Matrix{Int}, activeIndices::Vector{Int}, dualWeight::Vector{Oscar.TropicalSemiringElem{minOrMax}}, cellType::Symbol, ::typeof(max)=max)

Create a dual cell of the given type, using the max convention, with given ambient and active support.

# Arguments
- `ambientDualSupport::Matrix{Int}`: The ambient support of the dual cell, with columns as points.
- `activeIndices::Vector{Int}`: The active support of the dual cell, corresponding to indices of columns of the ambient support.
- `dualWeight::Vector{Oscar.TropicalSemiringElem{minOrMax}}`: The dual vector giving rise to this dual cell.
- `cellType::Symbol`: The type of the dual cell, must be one of :hypersurface, :linear, or :inverted_linear.
- `::typeof(max)`: The max convention.

# Throws
An error, as max tropical dual cells are currently unsupported.

"""
function dual_cell(ambientDualSupport::Matrix{Int}, activeIndices::Vector{Int}, dualWeight::Vector{Oscar.TropicalSemiringElem{typeof(max)}}, cellType::Symbol, ::typeof(max))
    error("max tropical dual cells currently unsupported")
end

"""
    check_dual_cell_inputs(ambientDualSupport::Matrix{Int}, activeIndices::Vector{Int}, cellType::Symbol, pluecker_vector::Union{TropicalPlueckerVector, Nothing}=nothing)

Check that the inputs to create a dual cell are valid.

# Arguments
- `ambientDualSupport::Matrix{Int}`: The ambient support of the dual cell.
- `activeIndices::Vector{Int}`: The active support of the dual cell, corresponding to indices of columns of the ambient support.
- `cellType::Symbol`: The type of the dual cell, must be one of :hypersurface, :linear, or :inverted_linear.
- `pluecker_vector::Union{TropicalPlueckerVector, Nothing}`: The pluecker vector of the dual cell, required for linear and inverted linear dual cells.

# Throws
An error if the inputs are invalid.
"""
function check_dual_cell_inputs(ambientDualSupport::Matrix{Int}, activeIndices::Vector{Int}, cellType::Symbol)
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
    return dim(convex_hull(active_support(m)))
end 

# handles printing of a dual cell
function Base.show(io::IO, m::DualCell)
    print(io, "Dual cell of type $(typeof(m).parameters[1]) supported on vertices $(m.activeIndices)")
end

function ambient_support(m::DualCell)
    return m.ambientDualSupport
end

function ambient_dim(m::DualCell)
    return size(ambient_support(m), 2)
end

function active_indices(m::DualCell)
    return m.activeIndices
end

function active_support(m::DualCell)
    return ambient_support(m)[active_indices(m)]
end

function dual_weight(m::DualCell)
    return m.dualWeight
end

function convention(m::DualCell{cellType, typeof(min)}) where cellType <: DualType
    return min
end

function convention(m::DualCell{cellType, typeof(max)}) where cellType <: DualType
    return max
end

function center(m::DualCell)
    return QQ.([sum(active_support(m), dims=1)...]) / length(active_indices(m))
end

import Oscar.convex_hull
function convex_hull(m::DualCell)
    return convex_hull(active_support(m))
end

function dual_facets(m::DualCell)
    cellHull = convex_hull(active_support(m))
    properDualFaces = DualCell[]
    for face in faces(cellHull, dim(cellHull) - 1)
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
function dual_cells(S::DualSupport{<:DualType}, c::Vector{<:Oscar.TropicalSemiringElem})

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