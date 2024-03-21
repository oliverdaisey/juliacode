include("../main.jl")

"""
This example demonstrates our framework for intersecting a tropical linear space and a tropical hypersurface.

We want to study how the intersection of the following tropical spaces behaves:
- (1) the parametric tropical hypersurface V(x1*x2*x3 + w*x0^3) in RR^4
- (2) the tropical linear space arising from the Pluecker vector in TT^(4c2) with all zeroes

The homotopy we follow initialises the system at w=-3 and moves to w=3, where inbetween (at w = 0) the number of intersection points changes.
"""

# set up the ambient rings
TT = tropical_semiring()
R, (x0, x1, x2, x3) = TT["x0", "x1", "x2", "x3"]

# write down the data for the tropical spaces
f_start = x1 * x2 * x3 + TT(-3)*x0^3
f_support = DualSupport{Hypersurface}(generate_support(f_start)) # helper function to generate the support from a polynomial
p_support = DualSupport{Linear}([1 1 0 0; 1 0 1 0; 1 0 0 1; 0 1 1 0; 0 1 0 1; 0 0 1 1]) # points of the matroid polytope projected down
# loopless facets: (1,2,3), (2,3,4), (1,3,4), (1,2,4), (3,5,6)

# write down the paths in the dual space
f_nodes = [TT.([0, -3]), TT.([0, 3])] # a path is given by a sequence of nodes
p_nodes = [TT.([0, 0, 0, 0, 0, 0])] # the tropical linear space s not moving
f_path = DualPath{Hypersurface, typeof(min)}(f_nodes, f_support)
p_path = DualPath{Linear, typeof(min)}(p_nodes, p_support)
h = mixed_path_in_series([f_path, p_path]) # this generates the path through the big dual space of their intersection

# create the starting dual cells and their mixed cell
f_start_dual = DualCell{Hypersurface, typeof(min)}(f_support, [1, 2], f_nodes[1])
p_start_dual = DualCell{Linear, typeof(min)}(p_support, [1, 2, 3], p_nodes[1]) # indices 1,2,3 correspond to a loopless facet
s = mixed_cell([f_start_dual, p_start_dual])

# create the mixed cell tracker to kick off our example
s_tracker = mixed_cell_tracker(h, s)

# verify that the mixed cell is correct
println("Tropical transverse intersection point = $(stable_intersection_point(s_tracker.mixed_cell))")

# verify the drift is in the correct direction
drift = tropical_drift(s_tracker)
println("Tropical drift = $(tropical_drift(s_tracker))")

# compute next breaking point
pt_of_interest, supports = next_point_of_interest(s_tracker)
println("next point of interest is $pt_of_interest")