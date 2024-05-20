include("../main.jl")

Kt, t = rational_function_field(QQ, "t")
R, z = Kt["z"]
nu = tropical_semiring_map(Kt, t)

Mtarget = matrix(R, [1 1 0 0 1 1; 0 0 1 1 1 1])
Mcurrent = matrix(R, [1 0 1 1 0 1; 0 1 0 1 1 1])
Mcurrent = matrix(R, rand(Int8, 2, 6))
Mtarget[:,3]=t.*Mtarget[:,3]
Mtarget[:,1]=t.*Mtarget[:,1]

M = Mtarget + z*Mcurrent

S = AbstractAlgebra.combinations(6,2)

entries = []
for s in S
    push!(entries, det(M[:,s]))
end

tropical_polynomial.(entries, Ref(nu))