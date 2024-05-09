include("../main.jl")

TT = tropical_semiring()
R, (x1, x2, x3, x4) = TT["x1", "x2", "x3", "x4"]

fStart = x1*x2 + 3*x3*x4
fTarget = x1*x2 + (-3)*x3*x4
fStartWeight = dual_weight(fStart)
fTargetWeight = dual_weight(fTarget)
fPath = dual_path([fStartWeight, fTargetWeight])

pWeight = dual_weight(Linear, [[-1, -1, 0, 0], [-1, 0, -1, 0], [-1, 0, 0, -1], [0, -1, -1, 0], [0 ,-1, 0, -1], [0, 0, -1, -1]], TT.([2, 0, 0, 0, 0, 0]))
pPath = dual_path(pWeight)

h = mixed_path_in_series([fPath, pPath])

# construct dual cells
fStartDual = dual_cell(fStartWeight)
pStartDual = dual_cell([[-1, 0, -1, 0], [-1, 0, 0, -1], [0, -1, -1, 0], [0 ,-1, 0, -1]], pWeight)

# construct mixed cell
s = mixed_cell([fStartDual, pStartDual])

# construct mixed cell tracker
t = mixed_cell_tracker(h, s)

# verify stable intersection point
stable_intersection_point(s)