include("../main.jl")

TT = tropical_semiring()
R, (x1, x2, x3) = TT["x1", "x2", "x3"]

f = x1*x2*x3 + 2

fSupport = DualSupport{Hypersurface}(f)
fNodes = [TT.([0, 2]), TT.([0, 6])] # a path is given by a sequence of nodes
fPath = dual_path(fNodes, fSupport)


pSupport = DualSupport{Linear}([1 1 0; 1 0 1; 1 0 0; 0 1 1; 0 1 0; 0 0 1])
pNodes = [TT.([2, 0, 0, 0, 0, 0])] # the tropical linear space is not moving
pPath = dual_path(pNodes, pSupport)

h = mixed_path_in_series([fPath, pPath])

fStartDual = DualCell{Hypersurface, typeof(min)}(fSupport, [1,2], fNodes[1])
pStartDual = DualCell{Linear, typeof(min)}(pSupport, [2, 3, 4, 5], pNodes[1])

s = mixed_cell([fStartDual, pStartDual])

t = mixed_cell_tracker(h, s)


