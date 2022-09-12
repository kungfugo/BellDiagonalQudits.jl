"""
    uniformBellSampler(n, d, object="magicSimplex", precision=10)

Create array of `n` uniformly distributed ``d^2`` Bell diagonal states represented as `CoordState` rounded to ``precision`` digits. 

Use `object="enclosurePolytope"` to create CoordStates in the enclosure polytope, having all ``coords \leq 1/d``.
"""
function uniformBellSampler(n, d, object="magicSimplex", precision=10)

    productParams = Array{CoordState}(undef, 0)
    nFound = 0

    while nFound < n

        u = zeros(d^2 + 1)
        u[2:(d^2)] = sort(vec(rand(Uniform(0, 1), 1, (d^2 - 1))))
        u[d^2+1] = 1

        coords = round.(diff(u), digits=precision)

        if object != "enclosurePolytope" || all(x -> x <= 1 / d, coords)
            push!(productParams, CoordState(coords, "UNKNOWN"))
            nFound += 1
        end

    end

    return productParams

end


"""
    createRandomCoordStates(nSamples, d, object="magicSimplex", precision=10, roundToSteps::Int=0, nTriesMax=10000000)

Return an array of `nSamples` `` d^2 `` dimensional CoordStates. 

Use the `object` to specify the coordinate ranges to [0,1] for 'magicSimplex' or [0, 1/d] for 'enclosurePolytope'. 
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

Return maximally entangled pure state of a bipartite system of dimension ``d^2``.
"""
function createBipartiteMaxEntangled(d)

    return proj(max_entangled(d * d))

end

"""
    weylOperator(d, k, l)

Return the ``(d,d)``- dimensional matrix representation of Weyl operator ``W_{k,l}``.
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

Return the tensor product ``W_{k,l} \\otimes \\mathbb{1}_d``.
"""
function getIntertwiner(d, k, l)
    w = weylOperator(d, k, l)
    Id = I(d)
    return w ⊗ Id
end

"""
    weylTrf(d, ρ, k, l)

Apply the ``(k,l)``-th Weyl transformation of dimension `d` to the density matrix `ρ`. Return ``W_{k,l} \\rho W_{k,l}^\\dagger``.
"""
function weylTrf(d, ρ, k, l)

    @assert size(ρ)[1] == d * d && size(ρ)[2] == d * d

    W = getIntertwiner(d, k, l)

    ρTrf = W * ρ * W'

    return ρTrf

end

"""
    createBasisStateOperators(d, bellStateOperator, precision)

Use maximally entangled Bell state `bellStateOperator` of dimension `d` to create Bell basis and return with corresponding flat and Weyl indices.
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

Return indexed Bell basis for `d` dimensions as `StandardBasis` rounded to `precision` digits.
"""
function createStandardIndexBasis(d, precision)::StandardBasis
    maxEntangled = createBipartiteMaxEntangled(d)
    standardBasis = StandardBasis(createBasisStateOperators(d, maxEntangled, precision))
    return standardBasis
end

"""
    createDensityState(coordState::CoordState, standardBasis::StandardBasis)

Return `DensityState`` containing the density matrix in computational basis based on `coordState` coordinates in Bell `standardBasis`.
"""
function createDensityState(coordState::CoordState, standardBasis::StandardBasis)::DensityState

    basisOps = map(x -> x[3], standardBasis.basis)
    densityState = DensityState(
        coordState.coords,
        Hermitian(genericVectorProduct(coordState.coords, basisOps)),
        coordState.eClass
    )

    return densityState

end

"""
    createWeylOperatorBasis(d)

Return vector of length ``d^2``, containing the Weyl operator basis for the ``(d,d)`` dimensionalmatrix space. 
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

Return vector of length ``d^4``, containing the product basis of two Weyl operator bases as basis for the ``(d^2,d^2)`` matrix space. 
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


"""
    createDimElementSubLattices(d)

Return all sublattices with `d` elements represented as vector of tuples in the ``d^2`` elements discrete phase space induced by Weyl operators.
"""
function createDimElementSubLattices(d)

    propDivs = getProperDivisors(d)

    allLattices = Array{Array{Tuple{Int,Int},1},1}(undef, 0)

    allPoints = Iterators.product(fill(0:(d-1), 2)...)

    for point in allPoints

        px = point[1]
        py = point[2]

        # Get horizontal line through point
        sSubLattice = Array{Tuple{Int,Int},1}(undef, 0)
        for s in (0:(d-1))
            sPoint = ((px + s) % d, py)
            push!(sSubLattice, sPoint)
        end
        push!(allLattices, sSubLattice)

        # Get vertical and diagonal lines through point
        for k in (0:(d-1))
            kSubLattice = Array{Tuple{Int,Int},1}(undef, 0)
            for s in (0:(d-1))
                ksPoint = ((px + k * s) % d, (py + s) % d)
                push!(kSubLattice, ksPoint)
            end
            push!(allLattices, kSubLattice)
        end

        # Get sublattices including point
        for b in propDivs
            for x in (0:b-1)
                xbSubLattice = Array{Tuple{Int,Int},1}(undef, 0)
                for u in (0:(d/b-1))
                    for v in (0:(b-1))
                        xuvbPoint = ((px + u * b + x * b) % d, (py + v * d / b) % d)
                        push!(xbSubLattice, xuvbPoint)
                    end
                end
                push!(allLattices, xbSubLattice)
            end
        end

    end

    # remove duplicates
    allLattices = unique(sort.(allLattices))
    return allLattices

end

"""
    createIndexSubLatticeState(standardBasis::StandardBasis, subLattice)

Return collection of `standardBasis` elements contributing to the state corresponding to the `sublattice`, coordinates in Bell basis and density matrix in computational basis.
"""
function createIndexSubLatticeState(standardBasis::StandardBasis, subLattice)

    lengthSub = length(subLattice)

    subBasisStates = Array{
        Tuple{
            Int,
            Array{Complex{Float64},2}
        },
        1}(undef, 0)

    #Go through all points and select for each the index and basis state
    for point in subLattice

        # find basis state that belongs to the given element of the sublattice
        foundIndex = findfirst(x -> x[2] == point, standardBasis.basis)
        push!(subBasisStates, (foundIndex, standardBasis.basis[foundIndex][3]))

    end

    indices = Array{Int,1}(undef, 0)
    subLatticeState = zeros(Complex{Float64}, size(standardBasis.basis[1][3]))

    # Collect indices and add basis state
    for state in subBasisStates
        push!(indices, state[1])
        subLatticeState = subLatticeState + state[2]
    end

    #Get corresponding coordinates
    coordinates = mapIndicesToNormCoord(indices, length(standardBasis.basis))

    # Collection of indices of basis elements which are equally mixed
    indexSubLatticeState = [indices, coordinates, 1 / lengthSub * subLatticeState]

    return indexSubLatticeState

end