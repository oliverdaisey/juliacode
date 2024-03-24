include("../main.jl")

"""
This example demonstrates our framework for intersecting an inverted tropical linear space and a tropical hypersurface.

We want to study how the intersection of the following tropical spaces behaves:
- (1) the parametric tropical hypersurface V(x1*x2*x3 + w*x0^3) in RR^4
- (2) the inverted tropical linear space arising from the Pluecker vector in TT^(4c2) with all zeroes

The homotopy we follow initialises the system at w=-3 and moves to w=3, where inbetween (at w = 0) the number of intersection points changes.
"""

# set up the ambient rings
TT = tropical_semiring()
R, (x0, x1, x2, x3) = TT["x0", "x1", "x2", "x3"]

# write down the data for the tropical spaces
fStart = x1 * x2 * x3 + TT(-3)*x0^3
fSupport = DualSupport{Hypersurface}(generate_support(fStart)) # helper function to generate the support from a polynomial
pSupport = DualSupport{InvertedLinear}([1 1 0 0; 1 0 1 0; 1 0 0 1; 0 1 1 0; 0 1 0 1; 0 0 1 1]) # points of the matroid polytope projected down
# loopless facets: (1,2,3), (2,3,4), (1,3,4), (1,2,4), (3,5,6)

# write down the paths in the dual space
fNodes = [TT.([0, -3]), TT.([0, 3])] # a path is given by a sequence of nodes
pNodes = [TT.([0, 0, 0, 0, 0, 0])] # the tropical linear space is not moving
fPath = DualPath{Hypersurface, typeof(min)}(fNodes, fSupport)
pPath = DualPath{InvertedLinear, typeof(min)}(pNodes, pSupport)
h = mixed_path_in_series([fPath, pPath]) # this generates the path through the big dual space of their intersection

# create the starting dual cells and their mixed cell
fStartDual = DualCell{Hypersurface, typeof(min)}(fSupport, [1, 2], fNodes[1])
pStartDual = DualCell{InvertedLinear, typeof(min)}(pSupport, [1, 2, 3], pNodes[1]) # indices 1,2,3 correspond to a loopless facet
s = mixed_cell([fStartDual, pStartDual])

# create the mixed cell tracker to kick off our example
sTracker = mixed_cell_tracker(h, s)

# verify that the mixed cell is correct
println("Tropical transverse intersection point = $(stable_intersection_point(sTracker.mixed_cell))")

# verify the drift is in the correct direction
drift = tropical_drift(sTracker)
println("Tropical drift = $(tropical_drift(sTracker))")

# compute all the mixed cells at the end
finalMixedCells = tropical_homotopy_continuation(sTracker)
println("Completed homotopy continuation, found $(length(finalMixedCells)) mixed cells")
for mixedCell in finalMixedCells
    println("Tropical transverse intersection point = $(stable_intersection_point(mixedCell))")
end