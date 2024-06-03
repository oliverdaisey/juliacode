include("../main.jl")

function starting_system(F::Vector{<:MPolyRingElem}, nu::TropicalSemiringMap; random_shift::Union{Vector{<:TropicalSemiringElem}, Nothing}=nothing)
    
    if length(F) == 1
        return starting_system_hypersurface(F, nu, random_shift=random_shift)
    elseif all(isequal(1),total_degree.(F))
        return starting_system_linear(F, nu, random_shift=random_shift)
    end

    error("Starting system not implemented for this case")
end

function starting_system_hypersurface(F::Vector{<:MPolyRingElem}, nu::TropicalSemiringMap; random_shift::Union{Vector{<:TropicalSemiringElem}, Nothing}=nothing)
    
    @assert length(F) == 1 "Only one polynomial is allowed"
    f = F[1]
    
    startingSystemSupport = sum(vars(f))+1
    if isnothing(random_shift)
        random_shift = tropical_semiring(nu).(rand(Int8, nvars(parent(f))))
    end
    randomLifts = random_lift.(Ref(nu),random_shift)
    h = evaluate(startingSystemSupport, randomLifts.*gens(parent(f)))

    d = total_degree(f)
    return [h^d], [h]

end

function starting_system_linear(F::Vector{<:MPolyRingElem}, nu::TropicalSemiringMap; random_shift::Union{Vector{<:TropicalSemiringElem}, Nothing}=nothing)
    
    @assert all(isequal(1),total_degree.(F)) "All polynomials must have total degree 1"

    R = parent(first(F))
    K = coefficient_ring(R)
    h = sum(gens(R)) + one(R)

    if isnothing(random_shift)
        random_shift = tropical_semiring(nu).(rand(Int8, nvars(R)))
    end

    randomLifts = gens(R).*random_lift.(Ref(nu),random_shift)

    startingSystem = elem_type(R)[]
    for i in 1:length(F)
        hRand = map_coefficients(c -> K(rand(Int8)), h)
        hRand = evaluate(hRand, randomLifts)
        push!(startingSystem, hRand)
    end

    return startingSystem, startingSystem

end

PolynomialSystem = Vector{<:MPolyRingElem}
StartingSystemData = Tuple{<:PolynomialSystem}
function starting_solution(startingSystems::Vector{StartingSystemData}, nu::TropicalSemiringMap) # add type for F

    system = vcat([it[2] for it in startingSystems]...)
    R = parent(first(system))

    @assert length(system) == ngens(R) "System is not square"
    A = zero_matrix(R, length(system), ngens(R) + 1)
    for (i,fi) in enumerate(system)
        for (j,xj) in enumerate(gens(R))
            A[i,j] = coeff(fi, xj)
        end
        A[i,ngens(R) + 1] = constant_coefficient(fi)
    end
    solution = nu.(kernel(A, side=:right))

    return [solution[i,1] / solution[end,1] for i in 1:(nrows(solution) - 1)]
end

function starting_data(partitionedSystem::Vector{PolynomialSystem}, nu::TropicalSemiringMap)

    startingSystems = starting_system.(partitionedSystem. Ref(nu))
    startingSolution = starting_solution(startingSystems, nu)

    for system in startingSystems
        if all(isequal(1), total_degree.(system[1]))
            flats = [findall(x -> x == val, solution) for val in unique(startingSolution)]
            if length(flats) < dim(ideal(system[1]))
                return starting_data(partitionedSystem, nu)
            end
        end
    end
end