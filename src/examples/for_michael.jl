include("../main.jl")

K, t = rational_function_field(QQ, "t")
nu = tropical_semiring_map(K, t)

R, (x1,x2,x3,x4,x5,x6) = polynomial_ring(K, ["x1","x2","x3","x4","x5","x6"])

f1 = x1 - x5*x6
f2 = x2 - x5*x6^2
f3 = x3 - x5*x6^3
f4 = x4 - x5*x6^4

F1 = [f1]
F2 = [f2]
F3 = [f3]
F4 = [f4]
F0 = [x1 + x2 + x5 + x6, x3 + x4 + x5 + x6]

partitionedSystem = [F1, F2, F3, F4, F0]
# activeIndices is not being computed correctly
startingSolution, startingSystems = starting_data(partitionedSystem, nu)

indicatorVectors = [
    [1, 1, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0],
    [1, 0, 0, 1, 0, 0],
    [1, 0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0, 1],
    [0, 1, 1, 0, 0, 0],
    [0, 1, 0, 1, 0, 0],
    [0, 1, 0, 0, 1, 0],
    [0, 1, 0, 0, 0, 1],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 0, 1, 0, 0, 1],
    [0, 0, 0, 1, 1, 0],
    [0, 0, 0, 1, 0, 1],
    [0, 0, 0, 0, 1, 1]
]

dotProducts = [sum(indicatorVector .* QQ.(startingSolution)) for indicatorVector in indicatorVectors]
# the minimum should be achieved multiple times