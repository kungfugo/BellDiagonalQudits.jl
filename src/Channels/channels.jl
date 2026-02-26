"""
    pauliweyl_channel(d, precision=10)

Construct a quantum channel in the Kraus representation using the generalized Bell basis.

The channel is defined by a set of d^2 Kraus operators which are the projectors 
onto the Bell basis states (Pauli-Weyl operators applied to a maximally entangled state). 

# Arguments
- `d::Int`: The dimension of the Hilbert space.
- `precision::Int`: Numerical precision for the internal operator construction.

# Returns
- `KrausOperators`: An object containing the collection of Bell basis projectors.
"""
function  pauliweyl_channel(d, precision=10)

    maxEntangled = create_bipartite_maxentangled(d)
    standardBasis = create_basis_state_operators(d, maxEntangled, precision)
    standardBellbasisOperators = map(x->x[3] , standardBasis)
    
    return KrausOperators(standardBellbasisOperators)

end

"""
    weyltwirl_channel_conj(d, precision=10)

Construct the bipartite Weyl-conjugate twirling channel in Kraus representation.

Generates \$d^2\$ Kraus operators of the form \$K_{k,l} = \\frac{1}{d} W_{k,l} \\otimes \\overline{W}_{k,l}\$, 
where \$W_{k,l}\$ are the generalized Pauli (Weyl) operators. This channel is typically 
used to symmetrize bipartite states into the isotropic state manifold.

# Arguments
- `d`: Subsystem dimension.
- `precision`: Rounding precision for the resulting matrices.
"""
function weyltwirl_channel_conj(d, precision=10)

    maxEntangled = create_bipartite_maxentangled(d)
    standardBasis = create_basis_state_operators(d, maxEntangled, precision)

    # Use same order of indices as used to create the standard basis.
    weylIndices = map(x->x[2], standardBasis)

    krausOps = Array{Complex{Float64},2}[]

    for (k,l) in weylIndices 
    
        W_kl = weyloperator(d,k,l)
        push!(krausOps, rounddigits(1/d * W_kl ⊗ conj(W_kl), precision))

    end

    return KrausOperators(krausOps)

end
 
"""
    weyltwirl_channel_trans(d, precision=10)

Construct the bipartite Weyl-transpose twirling channel in Kraus representation.

Generates \$d^2\$ Kraus operators of the form \$K_{k,l} = \\frac{1}{d} W_{k,l} \\otimes W_{k,l}^T\$.
While the conjugate twirl (\$W \\otimes \\overline{W}\$) targets isotropic states, this 
transpose twirl targets the space of Werner states, which are invariant under 
local unitary transformations of the form \$U \\otimes U\$.

# Arguments
- `d`: Subsystem dimension.
- `precision`: Rounding precision for the resulting matrices.
"""
function weyltwirl_channel_trans(d, precision=10)

    maxEntangled = create_bipartite_maxentangled(d)
    standardBasis = create_basis_state_operators(d, maxEntangled, precision)

    # Use same order of indices as used to create the standard basis.
    weylIndices = map(x->x[2], standardBasis)

    krausOps = Array{Complex{Float64},2}[]

    for (k,l) in weylIndices 
    
        W_kl = weyloperator(d,k,l)
        push!(krausOps, rounddigits(1/d * W_kl ⊗ transpose(W_kl), precision))

    end

    return KrausOperators(krausOps)

end

"""
    weyltwirl_channel_adj(d, precision=10)

Construct the bipartite Weyl-adjoint twirling channel in Kraus representation.

Generates \$d^2\$ Kraus operators of the form \$K_{k,l} = \\frac{1}{d} W_{k,l} \\otimes W_{k,l}^\\dagger\$.
This channel averages the input over bipartite Weyl unitaries where the second 
subsystem undergoes the adjoint transformation, often used for state 
symmetrization and identifying invariant subspaces.

# Arguments
- `d`: Subsystem dimension.
- `precision`: Rounding precision for the resulting matrices.
"""
function weyltwirl_channel_adj(d, precision=10)

    maxEntangled = create_bipartite_maxentangled(d)
    standardBasis = create_basis_state_operators(d, maxEntangled, precision)

    # Use same order of indices as used to create the standard basis.
    weylIndices = map(x->x[2], standardBasis)

    krausOps = Array{Complex{Float64},2}[]

    for (k,l) in weylIndices 
    
        W_kl = weyloperator(d,k,l)
        push!(krausOps, rounddigits(1/d * W_kl ⊗ W_kl', precision))

    end

    return KrausOperators(krausOps)

end

"""
    dephasing_channel(d)

Construct the completely dephasing channel in the computational basis.

This channel removes all off-diagonal elements (coherences) of a density matrix 
relative to the standard basis. The Kraus operators are the projectors 
\$P_i = |i\\rangle\\langle i|\$ for \$i \\in \\{1, \\dots, d\\}\$.

# Arguments
- `d::Int`: The dimension of the Hilbert space.

# Returns
- `KrausOperators`: A collection of \$d\$ rank-1 projection operators.
"""
function dephasing_channel(d)

    krausOps = Array{Complex{Float64},2}[]

    for i in 1:d
        K = proj(ket(i,d))
        push!(krausOps, K)
    end

    return(KrausOperators(krausOps))

end

"""
    local_dephasing_channel(dims, system)

Construct a quantum channel that applies complete dephasing to a specific subsystem.

For a multipartite system with dimensions `dims`, this function generates Kraus 
operators that decohere the subsystem at index `system`. The resulting Kraus 
operators are of the form I ⊗ ... ⊗ |i⟩⟨i| ⊗ ... ⊗ I.

# Arguments
- `dims`: A collection (tuple or vector) of dimensions for each subsystem.
- `system::Int`: The index of the subsystem to be dephased.

# Returns
- `KrausOperators`: A set of Kraus operators acting on the full Hilbert space.
"""
function local_dephasing_channel(dims, system)

    @assert system in 1:length(dims)

    krausOps = Array{Complex{Float64},2}[]

    dimSystem = dims[system]

    for i in 1:dimSystem 
        K = 1

        for ind in eachindex(dims)
            
            if ind == system 
                Ksys = proj(ket(i,dims[ind]))
            else
                Ksys = I(dims[ind])
            end

            K = K ⊗ Ksys

        end

        push!(krausOps,K)

    end

    return(KrausOperators(krausOps))

end

