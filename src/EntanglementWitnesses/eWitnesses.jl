"""
    createRandomWitnesses(standardBasis::StandardBasis, n)

Return array of `n` uniformly random `EntanglementWitness` (see `BellDiagonalQudits/structs.jl`) represented in `standardBasis`.
"""
function createRandomWitnesses(standardBasis::StandardBasis, n)::Array{EntanglementWitness}

    basisOps = map(x -> x[3], standardBasis.basis)
    D = length(basisOps)
    randomParams = rand(Uniform(-1, 1), n, D)

    wits = Array{EntanglementWitness}(undef, 0)
    for i in 1:n
        push!(wits, EntanglementWitness(randomParams[i, :], Hermitian(genericVectorProduct(randomParams[i, :], basisOps))))
    end

    return wits
end

"""
    createHalfSphericWitnesses(standardBasis::StandardBasis, n)

Return array of `n` uniformly random `EntanglementWitness` (see `BellDiagonalQudits/structs.jl`) on unit sphere represented in `standardBasis`.
"""
function createHalfSphericWitnesses(standardBasis::StandardBasis, n)::Array{EntanglementWitness}

    #Gaussian vectors with unit length in upper half plane are on unit sphere
    basisOps = map(x -> x[3], standardBasis.basis)
    D = length(basisOps)

    randomParams = Vector{Float64}[]

    for i in 1:n
        randomVector = normalize(rand(Normal(), D))
        if randomVector[1] < 0
            randomVector[1] = -randomVector[1]
        end

        push!(randomParams, randomVector)

    end

    wits = Array{EntanglementWitness}(undef, 0)
    for i in 1:n
        push!(wits, EntanglementWitness(randomParams[i], Hermitian(genericVectorProduct(randomParams[i], basisOps))))
    end

    return wits

end

"""
    getDirectFunctionsForWitTraceOptimization(wit::EntanglementWitness, d)

Return function and its negative that calculates trace of the given `wit` multiplied by a parameterized seperable state in `d` dimensions.
"""
function getDirectFunctionsForWitTraceOptimization(wit::EntanglementWitness, d)

    function getTrace(λ, witOpM, d)
        λ1 = λ[1:(2*(d-1))]
        λ2 = λ[(2*(d-1)+1):end]
        U1 = getCompRedParUnitaryFromVector(λ1, d)
        U2 = getCompRedParUnitaryFromVector(λ2, d)
        U = U1 ⊗ U2
        b = proj(ket(1, d^2))
        t = real(tr(witOpM * U * b * U'))
        return t
    end

    function opt_f(x)
        return getTrace(x, wit.operatorMatrix, d)
    end

    function opt_negf(x)
        return -getTrace(x, wit.operatorMatrix, d)
    end

    return (opt_f, opt_negf)

end


"""
    getWitnessExtrema(d, wit::EntanglementWitness, iterations, method, useConstrainedOpt=false)

Return optimization (see optimization.jl) results for lower and upper bound of `d` dimensional EntanglementWitness `wit` using `iterations` iterations and `method`.
"""
function getWitnessExtrema(
    d,
    wit::EntanglementWitness,
    iterations,
    method,
    useConstrainedOpt=false
)::Vector{Tuple{Float64,Vector{Float64}}}

    optFuncs = getDirectFunctionsForWitTraceOptimization(wit, d)
    optRes = directOptimization(
        optFuncs[1],
        optFuncs[2],
        method,
        d,
        iterations,
        useConstrainedOpt
    )

    return optRes

end

"""
    getBoundedEW(d, wit::EntanglementWitness, iterations, method=Optim.NelderMead, useConstrainedOpt=false)

Return BoundedEW in `d` dimensions based on EntanglementWitness `wit` and optimization of lower and upper bound.
"""
function getBoundedEW(
    d,
    wit::EntanglementWitness,
    iterations,
    method=Optim.NelderMead,
    useConstrainedOpt=false
)::BoundedEW

    witExtremizer = getWitnessExtrema(
        d,
        wit,
        iterations,
        method,
        useConstrainedOpt
    )

    b = proj(ket(1, d^2))
    totalMin = witExtremizer[1][1]
    totalMax = witExtremizer[2][1]

    minimizingParams = witExtremizer[1][2]
    minimizingParams1 = minimizingParams[1:(2*(d-1))]
    minimizingParams2 = minimizingParams[(2*(d-1)+1):end]
    minimizingU1 = getCompRedParUnitaryFromVector(minimizingParams1, d)
    minimizingU2 = getCompRedParUnitaryFromVector(minimizingParams2, d)
    minimizingU = minimizingU1 ⊗ minimizingU2
    minimizingDensityMatrix = Hermitian(minimizingU * b * minimizingU')

    maximizingParams = witExtremizer[2][2]
    maximizingParams1 = maximizingParams[1:(2*(d-1))]
    maximizingParams2 = maximizingParams[(2*(d-1)+1):end]
    maximizingU1 = getCompRedParUnitaryFromVector(maximizingParams1, d)
    maximizingU2 = getCompRedParUnitaryFromVector(maximizingParams2, d)
    maximizingU = maximizingU1 ⊗ maximizingU2
    maximizingDensityMatrix = Hermitian(maximizingU * b * maximizingU')

    return BoundedEW(
        wit.coords,
        totalMax,
        totalMin,
        maximizingDensityMatrix,
        minimizingDensityMatrix,
        iterations)

end
"""
    createRandomBoundedWits(
        d,
        standardBasis::StandardBasis,
        n,
        sphericalOnly::Bool,
        iterations::Integer,
        method=Optim.NelderMead,
        useConstrainedOpt=false
    )

Return array `n` BoundedEW with `standardBasis` coordinates uniformly distributed in [-1, 1] if `sphericalOnly` is false or uniformly distributed on sphere otherwise.

Use `iterations` runs to improve optimizatio with `method` impmlemented in Optim.jl.
"""
function createRandomBoundedWits(
    d,
    standardBasis::StandardBasis,
    n,
    sphericalOnly::Bool,
    iterations::Integer,
    method=Optim.NelderMead,
    useConstrainedOpt=false
)::Array{BoundedEW}

    if (length(standardBasis.basis) != d^2)
        throw("Dimension and basis do not match")
    end

    if (sphericalOnly)
        witnessCreator = createHalfSphericWitnesses
    else
        witnessCreator = createRandomWitnesses
    end

    randWits = witnessCreator(standardBasis, n)
    boundedEWs = map(x -> getBoundedEW(d, x, iterations, method, useConstrainedOpt), randWits)

    return boundedEWs

end

"""
    getBoundedCoordEw(bEw::BoundedEW)::BoundedCoordEW

Map BoundedEW  `bEw` to corresponding BoundedCoordEW.
"""
function getBoundedCoordEw(bEw::BoundedEW)::BoundedCoordEW

    return BoundedCoordEW(bEw.coords, bEw.upperBound, bEw.lowerBound, bEw.checkedIterations)

end