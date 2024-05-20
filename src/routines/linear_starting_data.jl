include("../main.jl")

function linear_starting_data(Mcurrent, Mtarget)
    
    @assert size(Mcurrent) == size(Mtarget) "Start and target matrix must have the same size"
    @assert base_ring(Mcurrent) == base_ring(Mtarget) "Start and target matrix must have the same base ring"
    @assert ngens(coefficient_ring(base_ring(Mcurrent))) == 1 "Base ring must be defined over a univariate rational function field"

    Kt = coefficient_ring(base_ring(Mcurrent))
    R, z = Kt["z"]
    t = gens(Kt)[1]
    nu = tropical_semiring_map(Kt, t)

    M = z*Mcurrent + Mtarget

    # for val(z) very negative, nu(M) should be constant modulo all ones vector

    S = AbstractAlgebra.combinations(ncols(M), nrows(M))

    entries = []
    for s in S
        push!(entries, det(M[:,s]))
    end

    tropicalisedSystem = tropical_polynomial.(entries, Ref(nu))

    # choose z really small and evaluate tropicalisedSystem
    smallZ = TT(-100000)
    
    dualVector = [entry(smallZ) for entry in tropicalisedSystem] 

    # normalise dualVector so that the first entry is 0
    dualVector = dualVector ./ dualVector[1]

    tHyp = tropical_hypersurface.(tropicalisedSystem)

    return dualVector, tropicalisedSystem

end

dualVector, tropicalisedSystem = linear_starting_data(Mcurrent, Mtarget)
println(tropicalisedSystem)
turningPoints = vertices.(polyhedral_complex.(tropical_hypersurface.(tropicalisedSystem)))

# points of tHyp are one-dimensional, so find the minimum
secondDualVector = TT.(minimum.([turningPoints[i][1] for i in 1:length(turningPoints)]))

