include("../main.jl")

TT = tropical_semiring()
R, (x1, x2, x3, x4) = TT["x1", "x2", "x3", "x4"]

f = x1*x2 + x3*x4
fSupport = DualSupport{Hypersurface}(f)
fNodes = [TT.([0,3]), TT.([0, 5])] # a path is given by a sequence of nodes
fPath = dual_path(fNodes, fSupport)


pSupport = DualSupport{Linear}([-1 -1 0 0; -1 0 -1 0; -1 0 0 -1; 0 -1 -1 0; 0 -1 0 -1; 0 0 -1 -1])
pNodes = [TT.([2, 0, 0, 0, 0, 0])] # the tropical linear space is not moving
pPath = dual_path(pNodes, pSupport)

h = mixed_path_in_series([fPath, pPath])

fStartDual = dual_cell(fSupport, [1,2], fNodes[1])
pStartDual = dual_cell(pSupport, [2, 3, 4, 5], pNodes[1])

s = mixed_cell([fStartDual, pStartDual])

t = mixed_cell_tracker(h, s)

# compute all the mixed cells at the end
finalMixedCells = tropical_homotopy_continuation(t)
println("Completed homotopy continuation, found $(length(finalMixedCells)) mixed cells")
for mixedCell in finalMixedCells
    println("Tropical transverse intersection point = $(stable_intersection_point(mixedCell))")
end