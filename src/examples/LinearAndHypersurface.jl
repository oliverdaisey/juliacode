include("../main.jl")

r"""
This example demonstrates our framework for intersecting a tropical linear space and a tropical hypersurface.

    We want to study how the intersection of the following tropical spaces behaves:
    - (1) the tropical hypersurface V(x1*x2*x3 + w*x0^3) in RR^4
    - (2) the tropical linear space arising from the Pluecker vector in TT^(4c2) with all zeroes

Since both spaces are invariant under translation by the all-ones vector, we project down to RR^3. This also gives a zero-dimensional intersection. 

The homotopy we follow initialises the system at w=3 and moves to w=-3, where inbetween (at w = 0) the number of intersection points changes.
"""

# set up the ambient rings
T = tropical_semiring()
R, (x1, x2, x3) = T["x1", "x2", "x3"]

# write down the data for the tropical spaces
f_start = x1*x2*x3 + T(3)
f_support = generate_support(f_start) # helper function to generate the support from a polynomial
p_support = [1 1 0; 1 0 1; 1 0 0; 0 1 1; 0 1 0; 0 0 1] # points of the matroid polytope projected down

# write down the paths in the dual space
f_nodes = [T.([0, 3]), T.([0, -3])] # a path is given by a sequence of nodes
p_nodes = [T.([0, 0, 0, 0, 0, 0])] # the tropical linear space is not moving
f_path = DualPath{DualPathHypersurface, typeof(min)}(f_nodes, f_support)
p_path = DualPath{DualPathLinear, typeof(min)}(p_nodes, p_support)
h = mixed_path_in_series([f_path, p_path]) # this generates the path through the big dual space of their intersection

# create the starting dual cells and their mixed cell
f_start_dual = DualCell{DualCellHypersurface, typeof(min)}(f_support, [1,2], f_nodes[1])
p_start_dual = DualCell{DualCellLinear, typeof(min)}(p_support, [1,2,3], p_nodes[1]) # indices 1,2,3 correspond to a loopless facet
s = mixed_cell([f_start_dual, p_start_dual])

# verify that the mixed cell is correct
println("Tropical transverse intersection point = $(stable_intersection_point(s, h, 1, QQFieldElem.(0)))")

# verify the drift is in the correct direction
println("Tropical drift = $(compute_drift(s, h, 1))")

# compute next breaking point
println("next breaking point at t = $(next_breaking_point(s, h, 1, QQFieldElem(0)))")