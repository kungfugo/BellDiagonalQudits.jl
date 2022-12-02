"""
    create_random_witnesses(standardBasis::StandardBasis, n)

Return array of `n` uniformly distributed random `EntanglementWitness` represented in Bell basis `standardBasis`.
"""
function create_random_witnesses(standardBasis::StandardBasis, n)::Array{EntanglementWitness}

    basisOps = map(x -> x[3], standardBasis.basis)
    D = length(basisOps)
    randomParams = rand(Uniform(-1, 1), n, D)

    wits = Array{EntanglementWitness}(undef, 0)
    for i in 1:n
        push!(wits, EntanglementWitness(randomParams[i, :], Hermitian(generic_vectorproduct(randomParams[i, :], basisOps))))
    end

    return wits
end

"""
    create_halfspheric_witnesses(standardBasis::StandardBasis, n)

Return array of `n` uniformly distributed random `EntanglementWitness` on unit sphere represented in Bell basis `standardBasis`.
"""
function create_halfspheric_witnesses(standardBasis::StandardBasis, n)::Array{EntanglementWitness}

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
        push!(wits, EntanglementWitness(randomParams[i], Hermitian(generic_vectorproduct(randomParams[i], basisOps))))
    end

    return wits

end

"""
    get_direct_functions_for_wit_traceoptimization(wit::EntanglementWitness, d)

Return the function and its negative that calculates ``tr \\rho`` `wit.coords`, the trace of the given witness `wit` multiplied by a parameterized seperable state ``\\rho`` in `d` dimensions.
"""
function get_direct_functions_for_wit_traceoptimization(wit::EntanglementWitness, d)

    function getTrace(λ, witOpM, d)
        λ1 = λ[1:(2*(d-1))]
        λ2 = λ[(2*(d-1)+1):end]
        U1 = get_comppar_unitary_from_parvector(λ1, d)
        U2 = get_comppar_unitary_from_parvector(λ2, d)
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
    get_witness_extrema(d, wit::EntanglementWitness, iterations, method, useConstrainedOpt=false)

Return optimization (see optimization.jl) results for lower and upper bound of `d` dimensional EntanglementWitness `wit` using `iterations` runs and Optim.jl optimization method `method`.
"""
function get_witness_extrema(
    d,
    wit::EntanglementWitness,
    iterations,
    method,
    useConstrainedOpt=false
)::Vector{Tuple{Float64,Vector{Float64}}}

    optFuncs = get_direct_functions_for_wit_traceoptimization(wit, d)
    optRes = direct_optimization(
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
    get_bounded_ew(d, wit::EntanglementWitness, iterations, method=Optim.NelderMead, useConstrainedOpt=false)

Return BoundedEW in `d` dimensions based on EntanglementWitness `wit` and `iterations` optimization runs of lower and upper bound for separable states.
"""
function get_bounded_ew(
    d,
    wit::EntanglementWitness,
    iterations,
    method=Optim.NelderMead,
    useConstrainedOpt=false
)::BoundedEW

    witExtremizer = get_witness_extrema(
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
    minimizingU1 = get_comppar_unitary_from_parvector(minimizingParams1, d)
    minimizingU2 = get_comppar_unitary_from_parvector(minimizingParams2, d)
    minimizingU = minimizingU1 ⊗ minimizingU2
    minimizingDensityMatrix = Hermitian(minimizingU * b * minimizingU')

    maximizingParams = witExtremizer[2][2]
    maximizingParams1 = maximizingParams[1:(2*(d-1))]
    maximizingParams2 = maximizingParams[(2*(d-1)+1):end]
    maximizingU1 = get_comppar_unitary_from_parvector(maximizingParams1, d)
    maximizingU2 = get_comppar_unitary_from_parvector(maximizingParams2, d)
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
    create_random_bounded_ews(
        d,
        standardBasis::StandardBasis,
        n,
        sphericalOnly::Bool,
        iterations::Integer,
        method=Optim.NelderMead,
        useConstrainedOpt=false
    )

Return array of `n` `BoundedEW` with ``d^2`` `standardBasis` coordinates uniformly distributed in [-1, 1] if `sphericalOnly` is false or uniformly distributed on unit sphere otherwise.

Use `iterations` runs to improve optimizatio with Optim.jl optimization method `method`.
"""
function create_random_bounded_ews(
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
        witnessCreator = create_halfspheric_witnesses
    else
        witnessCreator = create_random_witnesses
    end

    randWits = witnessCreator(standardBasis, n)
    boundedEWs = map(x -> get_bounded_ew(d, x, iterations, method, useConstrainedOpt), randWits)

    return boundedEWs

end

"""
    get_bounded_coordew(bEw::BoundedEW)::BoundedCoordEW

Map BoundedEW `bEw` to corresponding `BoundedCoordEW`.
"""
function get_bounded_coordew(bEw::BoundedEW)::BoundedCoordEW

    return BoundedCoordEW(bEw.coords, bEw.upperBound, bEw.lowerBound, bEw.checkedIterations)

end