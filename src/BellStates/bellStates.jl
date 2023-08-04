"""
    uniform_bell_sampler(n, d, object=:magicSimplex, precision=10)

Create array of `n` uniformly distributed ``d^2`` Bell diagonal states represented as `CoordState` rounded to `precision` digits. 

Use `object=:enclosurePolytope` to create CoordStates in the enclosure polytope, having all ``coords \\leq 1/d``.
"""
function uniform_bell_sampler(n, d, object=:magicSimplex, precision=10)

    productParams = Array{CoordState}(undef, 0)
    nFound = 0

    while nFound < n

        u = zeros(d^2 + 1)
        u[2:(d^2)] = sort(vec(rand(Uniform(0, 1), 1, (d^2 - 1))))
        u[d^2+1] = 1

        coords = round.(diff(u), digits=precision)

        if object != :enclosurePolytope || all(x -> x <= 1 / d, coords)
            push!(productParams, CoordState(coords, "UNKNOWN"))
            nFound += 1
        end

    end

    return productParams

end


"""
    create_random_coordstates(nSamples, d, object=:magicSimplex, precision=10, roundToSteps::Int=0, nTriesMax=10000000)

Return an array of `nSamples` `` d^2 `` dimensional CoordStates. 

Use the `object` to specify the coordinate ranges to [0,1] for 'magicSimplex' or [0, 1/d] for 'enclosurePolytope'. 
If `roundToSteps` > 0, round the coordinates to the vertices that divide the range in `roundToSteps`` equally sized sections.
Be aware that the resulting distribution of points is generally not uniform.
"""
function create_random_coordstates(nSamples, d, object=:magicSimplex, precision=10, roundToSteps::Int=0)::Array{CoordState}

    productParams = Array{CoordState}(undef, 0)
    nFound = 0

    if object != :enclosurePolytope

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
    create_bipartite_maxentangled(d)

Return maximally entangled pure state of a bipartite system of dimension ``d^2``.
"""
function create_bipartite_maxentangled(d)

    return proj(max_entangled(d * d))

end

"""
    weyloperator(d, k, l)

Return the ``(d,d)``- dimensional matrix representation of Weyl operator ``W_{k,l}``.
"""
function weyloperator(d, k, l)

    weylKetBra = zeros(Complex{Float64}, d, d)

    for j in (0:d-1)
        weylKetBra = weylKetBra + exp((2 / d * π * im) * j * k) * ketbra(j + 1, mod(j + l, d) + 1, d, d)
    end

    return weylKetBra

end

"""
    get_intertwiner(d, k, l)

Return the tensor product ``W_{k,l} \\otimes \\mathbb{1}_d``.
"""
function get_intertwiner(d, k, l)
    w = weyloperator(d, k, l)
    Id = I(d)
    return w ⊗ Id
end

"""
    weyltrf(d, ρ, k, l)

Apply the ``(k,l)``-th Weyl transformation of dimension `d` to the density matrix `ρ`. Return ``W_{k,l} \\rho W_{k,l}^\\dagger``.
"""
function weyltrf(d, ρ, k, l)

    @assert size(ρ)[1] == d * d && size(ρ)[2] == d * d

    W = get_intertwiner(d, k, l)

    ρTrf = W * ρ * W'

    return ρTrf

end

"""
    create_basis_state_operators(d, bellStateOperator, precision)

Use maximally entangled Bell state `bellStateOperator` of dimension `d` to create Bell basis and return with corresponding flat and Weyl indices.
"""
function create_basis_state_operators(d, bellStateOperator, precision)

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
        P_kl = Hermitian(round.(weyltrf(d, bellStateOperator, k, l), digits=precision))

        push!(basisStates, (index, value, P_kl))

    end
    return basisStates
end

"""
    create_standard_indexbasis(d, precision)

Return indexed Bell basis for `d` dimensions as `StandardBasis` rounded to `precision` digits.
"""
function create_standard_indexbasis(d, precision)::StandardBasis
    maxEntangled = create_bipartite_maxentangled(d)
    standardBasis = StandardBasis(create_basis_state_operators(d, maxEntangled, precision))
    return standardBasis
end

"""
    create_densitystate(coordState::CoordState, standardBasis::StandardBasis)

Return `DensityState`` containing the density matrix in computational basis based on `coordState` coordinates in Bell `standardBasis`.
"""
function create_densitystate(coordState::CoordState, standardBasis::StandardBasis)::DensityState

    basisOps = map(x -> x[3], standardBasis.basis)
    densityState = DensityState(
        coordState.coords,
        Hermitian(generic_vectorproduct(coordState.coords, basisOps)),
        coordState.eClass
    )

    return densityState

end

"""
    create_weyloperator_basis(d)

Return vector of length ``d^2``, containing the Weyl operator basis for the ``(d,d)`` dimensionalmatrix space. 
"""
function create_weyloperator_basis(d)::Vector{Array{Complex{Float64},2}}

    weylOperatorBasis = []
    for i in (0:(d-1))
        for j in (0:(d-1))
            push!(weylOperatorBasis, weyloperator(d, i, j))
        end
    end

    return weylOperatorBasis

end

"""
    create_bipartite_weyloperator_basis(d)

Return vector of length ``d^4``, containing the product basis of two Weyl operator bases as basis for the ``(d^2,d^2)`` matrix space. 
"""
function create_bipartite_weyloperator_basis(d)::Vector{Array{Complex{Float64},2}}

    weylOperatorBasis = create_weyloperator_basis(d)

    bipartiteWeylOperatorBasis = []

    for k in weylOperatorBasis
        for l in weylOperatorBasis
            push!(bipartiteWeylOperatorBasis, k ⊗ l)
        end
    end

    return bipartiteWeylOperatorBasis

end


"""
    create_dimelement_sublattices(d)

Return all sublattices with `d` elements represented as vector of tuples in the ``d^2`` elements discrete phase space induced by Weyl operators.
"""
function create_dimelement_sublattices(d)

    propDivs = get_properdivisors(d)

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
    create_index_sublattice_state(standardBasis::StandardBasis, subLattice)

Return collection of `standardBasis` elements contributing to the state corresponding to the `sublattice`, coordinates in Bell basis and density matrix in computational basis.
"""
function create_index_sublattice_state(standardBasis::StandardBasis, subLattice)

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
    coordinates = map_indices_to_normcoords(indices, length(standardBasis.basis))

    # Collection of indices of basis elements which are equally mixed
    indexSubLatticeState = [indices, coordinates, 1 / lengthSub * subLatticeState]

    return indexSubLatticeState

end

"""
    altbell_creatingoperator(d, k, l, α)

Return the ``(d,d)``- dimensional matrix representation of an alternative operator that is used to create Bell states ``V_{k,l}``.
``V_{k,l} ⊗ Id Ω_00`` creates ``1/√d * ∑ω^(ks)*α_(s,l) ket{s-l} ⊗ ket{s}``
"""
function altbell_creatingoperator(d, k, l, α::Matrix{Complex{Float64}})

    operatorKetBra = zeros(Complex{Float64}, d, d)

    @assert size(α) == (d, d)
    #Normalize
    α = map(x -> x / abs(x), α)


    for j in (0:d-1)
        operatorKetBra = operatorKetBra + (
            exp((2 / d * π * im) * (j + l) * k)
            * α[mod(j + l, d)+1, l+1]
            * ketbra(j + 1, mod(j + l, d) + 1, d, d)
        )
    end

    return operatorKetBra

end

"""
    create_altbellstate(state(d,k,l,α, returnDensity)

Return an alternative Bell state in the computational basis. 
Return denisty matrix unless returnDensity=false, in which case return state vector.
"""
function create_altbellstate(d, k, l, α::Matrix{Complex{Float64}}, returnDensity::Bool=true)

    V_kl = altbell_creatingoperator(d, k, l, α)
    IV_kl = V_kl ⊗ I(d)
    me = max_entangled(d * d)
    if returnDensity
        return (proj(IV_kl * me))
    else
        return (IV_kl * me)
    end
end

"""
    create_alt_indexbasis(d, l, α, precision)

Return alternative indexed Bell basis for `d` dimensions as `StandardBasis` rounded to `precision` digits.
Elements of secondary index `l` are replaced by alternative Bell state defined by `d`-element vector of phases `α`.
"""
function create_alt_indexbasis(d, α::Matrix{Complex{Float64}}, precision)::StandardBasis

    altIndexbasis = create_standard_indexbasis(d, precision)

    for (k, l) in Iterators.product(fill(0:(d-1), 2)...)

        #Weyl transformation to basis states
        P_alt_kl = rounddigits(Hermitian(create_altbellstate(d, k, l, α)), precision)

        #Find corresponding standard basis element
        elIndex = findfirst(x -> x[2] == (k, l), altIndexbasis.basis)

        #Replace element 
        altIndexbasis.basis[elIndex] = (elIndex, (k, l), P_alt_kl)

    end

    return (altIndexbasis)

end