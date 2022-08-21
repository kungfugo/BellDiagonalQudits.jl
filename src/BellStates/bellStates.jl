"""
    createRandomCoordStates(nSamples, d, object="magicSimplex", precision=10, roundToSteps::Int=0, nTriesMax=10000000)

Return an array of `nSamples` `d^2` dimensional CoordStates. 

Use the `object` to specify the coordinate ranges to [0,1] for 'magicSimplex' or [0, 1/d] for `enclosurePolytope`. 
If `roundToSteps` > 0, round the coordinates to the vertices that divide the range in `roundToSteps`` equally sized sections.
Be aware that the resulting distribution of points is generally not uniform.
"""
function createRandomCoordStates(nSamples, d, object="magicSimplex", precision=10, roundToSteps::Int=0)::Array{CoordState}

    productParams = Array{CoordState}(undef, 0)
    nFound = 0

    if object != "enclosurePolytope"

        while (nFound < nSamples)

            # only parameters which build a convex sum relevant
            c = vec(rand(Uniform(0, 1), 1, (d^2 - 1)))

            # If wanted, move random point to nearest step-point
            if roundToSteps > 0
                c = round.(roundToSteps * c) // roundToSteps
            end

            c = round.(c, digits=precision)

            if sum(c) <= 1
                push!(c, round(1 - sum(c), digits=precision))
                paramState = CoordState(c, "UNKNOWN")
                push!(productParams, paramState)
                nFound += 1
            end

        end

    else

        while (nFound < nSamples)

            c = vec(rand(Uniform(0, 1 / d), 1, (d^2 - 1)))

            if roundToSteps > 0
                c = round.(roundToSteps * d * c) // (roundToSteps * d)
            end

            c = round.(c, digits=precision)

            if (1 - 1 // d) <= sum(c) <= 1
                push!(c, round(1 - sum(c), digits=precision))
                paramState = CoordState(c, "UNKNOWN")
                push!(productParams, paramState)
                nFound += 1
            end

        end

    end

    return productParams

end

"""
    createBipartiteMaxEntangled(d)

Return maximally entangled pure state of a bipartite system of dimension `d`^2
"""
function createBipartiteMaxEntangled(d)

    maxE = proj(max_entangled(d * d))
    return maxE

end

"""
    weylOperator(d, k, l)

Return the (`d`,`d`)-dimensional matrix representation of Weyl operator with the indices `k` and `l`.
"""
function weylOperator(d, k, l)

    weylKetBra = zeros(Complex{Float64}, d, d)

    for j in (0:d-1)
        weylKetBra = weylKetBra + exp((2 / d * π * im) * j * k) * ketbra(j + 1, mod(j + l, d) + 1, d, d)
    end

    return weylKetBra

end

"""
    getIntertwiner(d, k, l)

Return the tensor product of the (`k`,`l`)-th Weyl operator and unity of dimension `d`.
"""
function getIntertwiner(d, k, l)
    w = weylOperator(d, k, l)
    Id = I(d)
    return w ⊗ Id
end

"""
    weylTrf(d, ρ, k, l)

Apply and return the (`k`,`l`)-th Weyl transformation of dimension `d` to the density matrix `ρ`.`
"""
function weylTrf(d, ρ, k, l)

    @assert size(ρ)[1] == d * d && size(ρ)[2] == d * d

    W = getIntertwiner(d, k, l)

    ρTrf = W * ρ * W'

    return ρTrf

end

"""
    createBasisStateOperators(d, bellStateOperator, precision)

Use `bellStateOperator` of dimension `d` to create Bell basis and return with corresponding flat and product indices.
"""
function createBasisStateOperators(d, bellStateOperator, precision)

    mIt = Iterators.product(fill(0:(d-1), 2)...)

    basisStates = Array{
        Tuple{
            Int,
            Tuple{Int,Int},
            Array{Complex{Float64},2}
        },
        1}(undef, 0)

    for (index, value) in enumerate(mIt)

        k = value[1]
        l = value[2]

        #Weyl transformation to basis states
        P_kl = Hermitian(round.(weylTrf(d, bellStateOperator, k, l), digits=precision))

        push!(basisStates, (index, value, P_kl))

    end
    return basisStates
end

"""
    createStandardIndexBasis(d, precision)

Return indexed Bell basis as type `StandardBasis` (see 'BellDiagonalQudits/src/structs.jl') with precision `precision` for `d` dimensions.
"""
function createStandardIndexBasis(d, precision)::StandardBasis
    maxEntangled = createBipartiteMaxEntangled(d)
    standardBasis = StandardBasis(createBasisStateOperators(d, maxEntangled, precision))
    return standardBasis
end

"""
    createDensityState(coordState::CoordState, indexBasis::StandardBasis)

Return DensityState (see 'BellDiagonalQudits/src/structs.jl') containing the density matrix in computational basis based on `coordState` coordinates in Bell `indexBais`.
"""
function createDensityState(coordState::CoordState, indexBasis::StandardBasis)::DensityState

    basisOps = map(x -> x[3], indexBasis.basis)
    densityState = DensityState(
        coordState.coords,
        Hermitian(genericVectorProduct(coordState.coords, basisOps)),
        coordState.eClass
    )

    return densityState

end

"""
    createWeylOperatorBasis(d)

Return vector of length `d`^2, containing the Weyl operator basis for the (`d`,`d`) matrix space. 
"""
function createWeylOperatorBasis(d)::Vector{Array{Complex{Float64},2}}

    weylOperatorBasis = []
    for i in (0:(d-1))
        for j in (0:(d-1))
            push!(weylOperatorBasis, weylOperator(d, i, j))
        end
    end

    return weylOperatorBasis

end

"""
    createBipartiteWeylOpereatorBasis(d)

Return vector of length `d`^4, containing the product basis of two Weyl operator bases as basis for the (`d`2,`d`2) matrix space. 
"""
function createBipartiteWeylOpereatorBasis(d)::Vector{Array{Complex{Float64},2}}

    weylOperatorBasis = createWeylOperatorBasis(d)

    bipartiteWeylOperatorBasis = []

    for k in weylOperatorBasis
        for l in weylOperatorBasis
            push!(bipartiteWeylOperatorBasis, k ⊗ l)
        end
    end

    return bipartiteWeylOperatorBasis

end