include("../main.jl")
include("../routines/generate_support.jl")

n=3
k=1
T = TropicalSemiring()
R, (x0, x1, x2, x3) = T["x0", "x1", "x2", "x3"]

# compute plueker indices
pluecker_indices = subsets(collect(1:(n+1)), k+1)
pluecker_vector = [0 for i in 1:length(pluecker_indices)]

f_start = x1*x2*x3 + T(3)*x0^3
f_target = x1*x2*x3 + T(-3)*x0^3

# create a DualCell for f
f_start_dual = DualCell{DualCellHypersurface, typeof(min)}(generate_support(f_start), [1,2])
# create the dual path
f_start_path = DualPath{DualPathHypersurface, typeof(min)}([[T(0),T(3)], [T(0), T(-3)]])

