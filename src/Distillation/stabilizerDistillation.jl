"""
    canonic_eigenbasis_weylprime(d, E)

Create array of eigenvectors-eigenvalue tuples for Weyl operator indexed by tuple `E` in dimension `d`.
"""
function canonic_eigenbasis_weylprime(d, E::Tuple{Int,Int})

    (a, b) = E
    ES = []

    if d == 2
        for n in 0:1
            if (a, b) == (0, 0)
                λ = 0
                ev_n = (modket(n, 2), λ)
            elseif b == 0
                λ = n
                ev_n = (modket(n * a, 2), λ)
            else
                λ = n + 1 / 2 * (mod(a * b, 2))
                ev_n = (1 / sqrt(2) * (modket(0, 2) + exp(2 / d * π * im*λ) * modket(1, 2)), λ)
            end
            push!(ES, ev_n)
        end
    else

        R, = residue_ring(ZZ, d)

        for λ in 0:(d-1)
            if (a, b) == (0, 0)
                ev_λ = (modket(λ, d), 0)
            elseif b == 0
                inv_a = Int(inv(R(a)).data)
                ev_λ = (modket(λ * inv_a, d), λ)
            else
                inv_b = Int(inv(R(b)).data)
                if d == 2
                    inv_2 = 0
                else
                    >
                    inv_2 = Int(inv(R(2)).data)
                end

                ev = 1 / sqrt(d) * modket(0, d)

                for j in 1:(d-1)
                    Γ_λj = Int(R(j * inv_b * λ - j * inv_b * (j * inv_b - 1) * inv_2 * a * b).data)
                    ev += 1 / sqrt(d) * exp(2 / d * π * im * Γ_λj) * modket(j, d)
                end
                ev_λ = (ev, λ)
            end
            push!(ES, ev_λ)
        end

    end

    return (ES)

end

"""
    create_canonic_enconding(g, d)

Return unitary matrix as encoding for generating stabilizer element `g` in demension`d`.
"""
function create_canonic_enconding(g::Vector{Tuple{Int,Int}}, d)

    @assert length(g) == 2
    (i_1, j_1) = g[1]
    (i_2, j_2) = g[2]

    ES1 = canonic_eigenbasis_weylprime(d, g[1])
    ES2 = canonic_eigenbasis_weylprime(d, g[2])

    if d == 2

        u = []
        for b in 0:1
            λ = mod(b + 1 / 2 * mod(i_1 * j_1 + i_2 * j_2, 2), 2)
            for k1 in 0:1
                (k1_v, k1_λ1) = (ES1[k1+1][1], ES1[k1+1][2])
                for k2 in 0:1
                    (k2_v, k2_λ2) = (ES2[k2+1][1], ES2[k2+1][2])

                    if (mod(k1_λ1 + k2_λ2, 2) == mod(λ, 2))
                        u_bk = k1_v ⊗ k2_v
                        push!(u, u_bk)
                    end

                end
            end
        end

    else

        u = []
        for b in 0:(d-1)
            for k1 in 0:(d-1)
                (k1_v, k1_b) = (ES1[k1+1][1], ES1[k1+1][2])
                for k2 in 0:(d-1)
                    (k2_v, k2_b) = (ES2[k2+1][1], ES2[k2+1][2])
                    if (mod(k1_b + k2_b, d) == b)
                        u_bk = k1_v ⊗ k2_v
                        push!(u, u_bk)
                    end
                end
            end
        end
    end


    @assert length(u) == d^2

    U = stack(u)
    return (U)

end

#Optional
function create_X_set_from_encodings(U, V)
    X_set = []

    for m in 0:d-1
        X_m = complex(zeros(d, d))
        for a in 0:(d-1), b in 0:(d-1)
            ma_index = findfirst(modket(m, d) ⊗ modket(a, d) .== 1)
            mb_index = findfirst(modket(m, d) ⊗ modket(b, d) .== 1)
            u_ma = U[:, ma_index]
            v_mb = V[:, mb_index]

            X_m += u_ma' * v_mb * modket(a, d) * modbra(b, d)
        end
        push!(X_set, X_m)
    end

    return (X_set)

end

"""
    get_s_for_generator(g, E)

Return s-value of error element `E`` for stabilzer generator `g` in dimension `d`.
"""
function get_s_for_generator(g::Vector{Tuple{Int,Int}}, E::Vector{Tuple{Int,Int}}, d)

    @assert length(g) == length(E)

    s = 0
    for n in eachindex(g)
        s += g[n][2] * E[n][1] - g[n][1] * E[n][2]
    end

    return (mod(s, d))

end

"""
    get_erroroperator_from_errorelement(E, d)

Return error operator in computational basis from error element `E` in dimension `d`.
"""
function get_erroroperator_from_errorelement(E::Vector{Tuple{Int,Int}}, d)

    e = first(E)
    w = weyloperator(d, e[1], e[2])

    for el in E[2:lastindex(E)]
        w = w ⊗ weyloperator(d, el[1], el[2])
    end

    return (w)

end


"""
    get_erroraction_from_encoding(U, W, bVec, sVec, d )

Create error action operator for an error `W` in encoding `U` for measurment outcomes `bVec` and `sVec` in dimension `d`.
"""
function get_erroraction_from_encoding(U, W, bVec, sVec, d)

    @assert size(U) == size(W)
    @assert length(bVec) == length(sVec) == 1
    b = bVec[1]
    s = sVec[1]

    Y = U' * W * U

    T = complex(zeros(d, d))
    for k in 0:(d-1), l in 0:(d-1)

        Tkl = (modbra(b + s, d) ⊗ modbra(k, d)) * Y * (modket(b, d) ⊗ modket(l, d))

        T += Tkl * (modket(k, d) ⊗ modbra(l, d))

    end

    return T
end

"""
    get_s_dict_for_generators(gVec, n, d)

Create Dict containing the s-values for each error element (keys) for stabilizer generating elements `gVec` for `n`-copies in dimension `d`.
"""
function get_s_dict_for_generators(gVec::Vector{Vector{Tuple{Int,Int}}}, n, d)

    @assert n == 2
    @assert length(gVec) == 1

    dict = Dict()

    for i1 in 0:(d-1), j1 in 0:(d-1), i2 in 0:(d-1), j2 in 0:(d-1)
        E = [(i1, j1), (i2, j2)]
        sVec = Int[]

        for g in gVec
            s = get_s_for_generator(g, E, d)
            push!(sVec, s)
        end

        dict[E] = sVec
    end

    return (dict)

end

"""
    get_codespace_errorops_for_meas(U, gVec, aVec, bVec)

For given encoding `U` of a stabilizer with generating elements `gVec` and measurement outcomes `aVec` and `bVec`, return Dict of error elements (keys) and error action operators (values).
"""
function get_codespace_errorops_for_meas(U, gVec::Vector{Vector{Tuple{Int,Int}}}, aVec::Vector{Int}, bVec::Vector{Int}, n, d)
    @assert length(aVec) == length(bVec) == 1
    @assert n == 2

    sVec = map(x -> mod(x, d), aVec - bVec)

    s_dict = get_s_dict_for_generators(gVec, n, d)
    codespace_errorops_dict = filter(x -> x[2] == sVec, s_dict)

    for errorkey in keys(codespace_errorops_dict)

        W_errorop = get_erroroperator_from_errorelement(errorkey, d)
        cs_errorop = get_erroraction_from_encoding(U, W_errorop, bVec, sVec, d)
        codespace_errorops_dict[errorkey] = cs_errorop

    end

    return (codespace_errorops_dict)
end

"""
    get_stdbasis_probabilities(ρ, stdBasis, precision)

Return probability of measuring a Bell `stdBasis` state for `ρ` in `precision`.
"""
function get_stdbasis_probabilities(ρ, stdBasis, precision=10)
    basis = stdBasis.basis

    probs = map(x -> (x[1], x[2], (tr(ρ * x[3]) |> real |> (y -> rounddigits(y, precision)))), basis)

    return (probs)

end

"""
    get_prob_for_E(ρ, E, stdBasis)

Return probability of measurment of `stdBasis` corresponding to error element `E` in state `ρ`.
"""
function get_prob_for_E(ρ, E::Vector{Tuple{Int,Int}}, stdbasis)
    stdbasisProbs = get_stdbasis_probabilities(ρ, stdbasis, 10)
    pE = 1
    for Em in E
        pEm = real(first(filter(x -> x[2] == Em, stdbasisProbs))[3])
        pE *= pEm
    end
    return (pE)
end

"""
    get_s_errorkeys(gVec, sVec, n, d)

Return all errors with `sValue` for `n`-cope stabilizer generators `gVec` in `d` dimenions.
"""
function get_s_errorkeys(gVec, sVec, n, d)
    @assert n == 2

    sDict = get_s_dict_for_generators(gVec, n, d)
    ES = keys(filter(x -> x[2] == sVec, sDict))

    return (ES)

end

"""
    get_prob_for_errorkeys(ρ, errorKeys, n, stdBasis)

Get total error probability for `n`-copy `errorKeys` defined by `stdBasis` in `rho`.
"""
function get_prob_for_errorkeys(ρ, errorKeys, n, stdbasis)

    @assert n == 2

    pS = 0
    for E in errorKeys
        pE = get_prob_for_E(ρ, E, stdbasis)
        pS += pE
    end

    return (pS)

end

"""
    get_abelian_subgroup_generators(d, n)

Return all generators of abelian subgroups for the `n`-copy weyl-Heisenberg group in `d` dimensions.
"""
function get_abelian_subgroup_generators(d, n)

    @assert n == 2
    @assert d ∈ (2, 3)

    potentialGenerators = Vector{Tuple{Int,Int}}[]
    generators = Vector{Tuple{Int,Int}}[]

    for i1 in 0:(d-1), j1 in 0:(d-1), i2 in 0:(d-1), j2 in 0:(d-1)

        if (i1, j1, i2, j2) !== (0, 0, 0, 0)
            E = [(i1, j1), (i2, j2)]
            push!(potentialGenerators, E)
        end
    end

    while !isempty(potentialGenerators)

        g = first(potentialGenerators)
        push!(generators, g)

        gg = add_errorkeys(g, g, d, n)

        filter!(x -> x != g, potentialGenerators)
        filter!(x -> x != gg, potentialGenerators)

    end

    return (generators)

end

"""
    get_coset_partition(gVec, n, d)

Return coset partition of error elements defined by stabilizer generating elements `gVec` of `n`-copies in `d` dimensions.
"""
function get_coset_partition(gVec, n, d)

    @assert n == 2
    @assert d ∈ (2, 3)
    @assert length(gVec) == 1
    @assert gVec[1] !== [(0, 0), (0, 0)]

    g = gVec[1]
    gg = add_errorkeys(g, g, d, n)

    errorKeys = Vector{Tuple{Int,Int}}[]


    for i1 in 0:(d-1), j1 in 0:(d-1), i2 in 0:(d-1), j2 in 0:(d-1)
        if (i1, j1, i2, j2) !== (0, 0, 0, 0)
            E = [(i1, j1), (i2, j2)]
            push!(errorKeys, E)
        end
    end

    cosets = Vector{Vector{Tuple{Int,Int}}}[]

    while !isempty(errorKeys)

        coset = Vector{Tuple{Int,Int}}[]

        E = first(errorKeys)
        push!(coset, E)

        for x in 1:(d-1)
            E = add_errorkeys(E, g, d, n)
            push!(coset, E)
        end

        push!(cosets, coset)

        filter!(x -> !(x ∈ coset), errorKeys)

    end

    return (cosets)

end

"""
    get_cosets_prob_s(ρ, gVec, n, d, stdbasis)

Return vector of tuples containing the cosets defined by stabilizer generating elements `gVec` for `n`-copies in`d` dimensions together with its s-value and probabilies in the `stdbasis`.
"""
function get_cosets_prob_s(ρ, gVec, n, d, stdbasis)

    s_dict = get_s_dict_for_generators(gVec, n, d)
    cosets = get_coset_partition(gVec, n, d)

    cosetsProbS = map(c -> (c, s_dict[c[1]], get_prob_for_errorkeys(ρ, c, n, stdbasis)), cosets)

    return (cosetsProbS)

end

"""
    get_s_prob_dict(ρ, gVec, n, d, stdbasis)

Return Dict of `stdbasis` probabilities (values) for each s-value (keys) in `n`-copy state `rho` of subsystem dimension `d`.
"""
function get_s_prob_dict(ρ, gVec, n, d, stdbasis)

    s_dict = get_s_dict_for_generators(gVec, n, d)
    s_values = unique(values(s_dict))

    s_prob_dict = Dict()

    for sVec in s_values
        ES = get_s_errorkeys(gVec, sVec, n, d)
        PS = get_prob_for_errorkeys(ρ, ES, n, stdbasis)
        s_prob_dict[sVec] = PS
    end

    return (s_prob_dict)
end

"""
    get_s_normalized_coset_probs(ρ, gVec, n, d, stdbasis)

Return array of tuples containing the cosets and s-values defined by generating `n`-copy stabilizer elements `gVec` 
and the coset probabilitiies, normalized coset probabilites and the s-value probabilies for the `stdbasis` in `d` dimensions.
"""
function get_s_normalized_coset_probs(ρ, gVec, n, d, stdbasis)

    cosetProbS = get_cosets_prob_s(ρ, gVec, n, d, stdbasis)
    sProbDict = get_s_prob_dict(ρ, gVec, n, d, stdbasis)
    #remove zero probability sVecs
    filter!(x -> x[2] > 0, sProbDict)
    filter!(x -> x[2] ∈ keys(sProbDict), cosetProbS)

    normCosetProbs = map(x -> (x[1], x[2], x[3], x[3] / sProbDict[x[2]], sProbDict[x[2]]), cosetProbS)


    return (sort(normCosetProbs, by=(x -> -x[4])))

end

"""
    get_max_fidelity_stabilizer_coset(ρ, n, d, stdbasis)

Returns stabilizer generator, coset of maximal normalized probability and corresponding s-value, coset probabily and s-value probabily for `n`-copy input state `ρ` in `d` dimensions with respect to `stdbasis`.
"""
function get_max_fidelity_stabilizer_coset(ρ, n, d, stdbasis)

    allGens = get_abelian_subgroup_generators(d, n)

    maxFidelityStabCosets = map(
        g -> (g, first(get_s_normalized_coset_probs(ρ, [g], n, d, stdbasis))),
        allGens
    )

    opt = sort(maxFidelityStabCosets, by=(x -> -x[2][4])) |> first


    return (opt[1], opt[2][1], opt[2][2], opt[2][4], opt[2][5])

end

"""
    stabilizer_routine((ρ, gVec, aVec, bVec, Uenc, n, d, stdbasis)

Return output state and probability of success for one iteration of the stabilizer distillation routine.

The routine is defined by the input state `ρ`, the stabilizer generating element `gVec`, the measurement outcomes `aVec` and `bVec`, 
the encoding operator `Uenc`, the number of copies `n`, the dimension of the subsystems `d` and the standard Bell bases `stdbasis`.
"""
function stabilizer_routine(ρ, gVec, aVec, bVec, Uenc, n, d, stdbasis::StandardBasis)

    @assert length(stdbasis.basis) == d^2
    @assert size(ρ) == (d^2, d^2)
    @assert n == 2
    @assert length(gVec) == length(aVec) == length(bVec) == 1


    sVec = aVec - bVec
    ES = get_s_errorkeys(gVec, sVec, n, d)
    PS = get_prob_for_errorkeys(ρ, ES, n, stdbasis)

    if rounddigits(PS, 10) == 0
        return (ρ)
    end

    TES = get_codespace_errorops_for_meas(Uenc, gVec, aVec, bVec, n, d)

    P00 = stdbasis.basis[1][3]

    ρ_out = complex(zeros(d^2, d^2))
    for E in ES
        pE = get_prob_for_E(ρ, E, stdbasis)
        TE = TES[E]
        ρ_out += pE * (
            (TE ⊗ I(d)) * P00 * (TE ⊗ I(d))'
        )
    end
    ρ_out = 1 / PS * ρ_out #Normalization
    return (ρ_out, PS)

end

"""
    FIMAX_routine(ρ, n, d, stdbasis)

Applies one iteration of FIMAX routine to the `n`-copy input state `'rho` in `d` dimensions. Returns the output state and the probabilies of success with respect to the `stdbasis`.
"""
function FIMAX_routine(ρ, n, d, stdbasis::StandardBasis)
    @assert n == 2
    @assert d ∈ (2, 3)

    maxFidelityStabilizerCoset = get_max_fidelity_stabilizer_coset(ρ, n, d, stdbasis)

    (g, cosetRep, sVec, successProb) = (x -> (x[1], x[2][1], x[3], x[5]))(maxFidelityStabilizerCoset)

    aVec = sVec
    bVec = aVec - sVec

    U_canonic = create_canonic_enconding(g, d)
    cosetRepOp = get_erroroperator_from_errorelement(cosetRep, d)
    cosetRepCodespaceErrorOp = get_erroraction_from_encoding(U_canonic, cosetRepOp, bVec, sVec, d)

    ρ_out_raw = stabilizer_routine(ρ, [g], aVec, bVec, U_canonic, n, d, stdbasis)[1]

    #Correction operation
    correctionOp = cosetRepCodespaceErrorOp'

    ρ_out = (correctionOp ⊗ I(d)) * ρ_out_raw * (correctionOp ⊗ I(d))'

    return (ρ_out, real(successProb))

end

"""
    P1_routine(ρ, n, d, stdbasis)

Applies one iteration of P1 routine to the `n`-copy input state `'rho` in `d` dimensions. Returns the output state and the probabilies of success with respect to the `stdbasis`.
"""
function P1_routine(ρ, n, d, stdbasis::StandardBasis)

    g = [(1, 0), (mod(-1, d), 0)]

    gVec = [g]
    aVec = [0]
    bVec = [0]

    U = create_canonic_enconding(g, d)
    (ρ, succProb) = stabilizer_routine(ρ, gVec, aVec, bVec, U, n, d, stdbasis)

    return (ρ, real(succProb))

end

"""
    P2_routine(ρ, n, d, stdbasis)

Applies one iteration of P2 routine to the `n`-copy input state `'rho` in `d` dimensions. Returns the output state and the probabilies of success with respect to the `stdbasis`.
"""
function P2_routine(ρ, n, d, stdbasis::StandardBasis)

    g = [(1, 0), (mod(-1, d), 0)]

    gVec = [g]
    aVec = [0]
    bVec = [0]

    U = create_canonic_enconding(g, d)

    bFT = (FT_OP(d) ⊗ conj(FT_OP(d)))

    ρ = bFT * ρ * bFT'
    (ρ, succProb) = stabilizer_routine(ρ, gVec, aVec, bVec, U, n, d, stdbasis)
    ρ = bFT * ρ * bFT'

    return (ρ, real(succProb))

end

"""
    P1_P2_routine(ρ, n, d, stdbasis)

Applies one iteration of P1_P2 routine to the `n`-copy input state `'rho` in `d` dimensions. Returns the output state and the probabilies of success with respect to the `stdbasis`.
"""
function P1_P2_routine(ρ, n, d, stdbasis::StandardBasis)

    basisDict = create_dictionary_from_basis(stdbasis)
    bellCoords = belldiagonal_projection(stdbasis, ρ).coords

    Z_indices = map(i -> basisDict[1][i], [(k, 0) for k in 0:(d-1)])
    X_indices = map(i -> basisDict[1][i], [(0, k) for k in 0:(d-1)])

    Z_error_sum = sum(bellCoords[Z_indices])
    X_error_sum = sum(bellCoords[X_indices])

    if Z_error_sum > X_error_sum
        (ρ, succProb) = P2_routine(ρ, n, d, stdbasis)
    else
        (ρ, succProb) = P1_routine(ρ, n, d, stdbasis)
    end

    return (ρ, real(succProb))
end

"""
    BBPSSW_routine(ρ, n, d, stdbasis)

Applies one iteration of BBPSSW routine to the `n`-copy input state `'rho` in `d` dimensions. Returns the output state and the probabilies of success with respect to the `stdbasis`.
"""
function BBPSSW_routine(ρ, n, d, stdbasis::StandardBasis)
    @assert n == 2

    basisDict = create_dictionary_from_basis(stdbasis)
    densityState = belldiagonal_projection(stdbasis, ρ)
    coords = densityState.coords

    depol_coords = depolarize_coords(coords)

    bellCoords(k, l) = depol_coords[basisDict[1][(k, l)]]
    newCoords = zeros(d^2)

    for x in 0:(d-1), y in 0:(d-1)

        i = basisDict[1][(x, y)]
        a_xy = 0.0

        for k in 0:(d-1)
            a_xy += bellCoords(k, y) * bellCoords(mod(x - k, d), y)
        end

        newCoords[i] = a_xy

    end

    succProb = sum(newCoords)

    if succProb == 0.0
        throw("successProb is zero")
    end

    normNowCoords = newCoords / succProb
    ρ_new = create_densitystate(CoordState(normNowCoords, "UNKNOWN"), stdbasis).densityMatrix

    return (ρ_new, real(succProb))

end

"""
    DEJMPS_routine(ρ, n, d, stdbasis)

Applies one iteration of DEJMPS routine to the `n`-copy input state `'rho` in `d` dimensions. Returns the output state and the probabilies of success with respect to the `stdbasis`.
"""
function DEJMPS_routine(ρ, n, d, stdbasis::StandardBasis)

    g = [(1, 0), (d - 1, 0)]

    gVec = [g]
    aVec = [0]
    bVec = [0]

    U = create_canonic_enconding(g, d)

    bFT = (FT_OP(d) ⊗ conj(FT_OP(d)))

    (ρ, succProb) = stabilizer_routine(ρ, gVec, aVec, bVec, U, n, d, stdbasis)
    ρ = bFT * ρ * bFT'

    return (ρ, real(succProb))

end

"""
    iterative_specific_stabilizer_protocol(ρ, targetFid, gVec, aVec, bVec, U, n, d, stdBasis, maxIts)

Applies iterations of the stabilizer routine to the `n`-copy input state `'rho` in `d` dimensions for stabilizer generating elements `gVec`, encoding `U` and measurement outcomes `aVec` and`bVec`.
Iterates until `targetFid` or maximal number of iterations `maxIts` is reached. 
Returns `distillable=true` if `targetFid` could be reached and fidelities and success probabilies with respect to the `stdbasis` for each iteration.
"""
function iterative_specific_stabilizer_protocol(ρ, targetFid, gVec, aVec, bVec, U, n, d, stdbasis::StandardBasis, maxIts=100)
    P_00 = stdbasis.basis[1][3]
    isPpt = isppt(ρ, d, 10)
    F_00 = fidelity(ρ, P_00)

    i = 0
    iterations = [i]
    fidelities = [F_00]
    successProbs = [0.0]


    while (F_00 <= targetFid) && !isPpt && (i < maxIts)

        i += 1

        (ρ, succProb) = stabilizer_routine(ρ, gVec, aVec, bVec, U, n, d, stdbasis)
        isPpt = isppt(ρ, d, 10)
        F_00 = fidelity(ρ, P_00)

        push!(iterations, i)
        push!(fidelities, F_00)
        push!(successProbs, succProb)

    end

    distillable = false
    if F_00 > targetFid && i < maxIts
        distillable = true
    elseif isPpt || i >= maxIts
        distillable = false
    end

    return (distillable, iterations, fidelities, successProbs)

end

"""
    iterative_FIMAX_protocol(ρ, targetFid, n, d, stdBasis, maxIts)

Applies iterations of the FIMAX routine to the `n`-copy input state `ρ` in `d` dimensions. 
Iterates until `targetFid` or maximal number of iterations `maxIts` is reached. 
Returns `distillable=true` if `targetFid` could be reached and fidelities and success probabilies with respect to the `stdbasis` for each iteration.
"""
function iterative_FIMAX_protocol(ρ, targetFid, n, d, stdbasis, maxIts=100)
    P_00 = stdbasis.basis[1][3]
    isPpt = isppt(ρ, d, 10)
    F_00 = fidelity(ρ, P_00)

    i = 0
    iterations = [i]
    fidelities = [F_00]
    successProbs = [0.0]


    while (F_00 <= targetFid) && !isPpt && (i < maxIts)

        i += 1

        (ρ, succProb) = FIMAX_routine(ρ, n, d, stdbasis)
        isPpt = isppt(ρ, d, 10)
        F_00 = fidelity(ρ, P_00)

        push!(iterations, i)
        push!(fidelities, F_00)
        push!(successProbs, succProb)

    end

    distillable = false
    if F_00 > targetFid && i < maxIts
        distillable = true
    elseif isPpt || i >= maxIts
        distillable = false
    end

    return (distillable, iterations, fidelities, successProbs)

end

"""
    iterative_P1_P2_protocol(ρ, targetFid, n, d, stdBasis, maxIts)

Applies iterations of the P1_P2 routine to the `n`-copy input state `ρ` in `d` dimensions. 
Iterates until `targetFid` or maximal number of iterations `maxIts` is reached. 
Returns `distillable=true` if `targetFid` could be reached and fidelities and success probabilies with respect to the `stdbasis` for each iteration.
"""
function iterative_P1_P2_protocol(ρ, targetFid, n, d, stdbasis, maxIts=100)
    P_00 = stdbasis.basis[1][3]
    isPpt = isppt(ρ, d, 10)
    F_00 = fidelity(ρ, P_00)

    i = 0
    iterations = [i]
    fidelities = [F_00]
    successProbs = [0.0]

    while (F_00 <= targetFid) && !isPpt && (i < maxIts)

        i += 1

        (ρ, succProb) = P1_P2_routine(ρ, n, d, stdbasis)
        isPpt = isppt(ρ, d, 10)
        F_00 = fidelity(ρ, P_00)

        push!(iterations, i)
        push!(fidelities, F_00)
        push!(successProbs, succProb)

    end

    distillable = false
    if F_00 > targetFid && i < maxIts
        distillable = true
    elseif isPpt || i >= maxIts
        distillable = false
    end

    return (distillable, iterations, fidelities, successProbs)

end

"""
    iterative_BBPSSW_protocol(ρ, targetFid, n, d, stdBasis, maxIts)

Applies iterations of the BBPSSW routine to the `n`-copy input state `ρ` in `d` dimensions. 
Iterates until `targetFid` or maximal number of iterations `maxIts` is reached. 
Returns `distillable=true` if `targetFid` could be reached and fidelities and success probabilies with respect to the `stdbasis` for each iteration.
"""
function iterative_BBPSSW_protocol(ρ, targetFid, n, d, stdbasis, maxIts=100)
    P_00 = stdbasis.basis[1][3]
    isPpt = isppt(ρ, d, 10)
    F_00 = fidelity(ρ, P_00)

    i = 0
    iterations = [i]
    fidelities = [F_00]
    successProbs = [0.0]

    while (F_00 <= targetFid) && !isPpt && (i < maxIts)

        i += 1

        (ρ, succProb) = BBPSSW_routine(ρ, n, d, stdbasis)
        isPpt = isppt(ρ, d, 10)
        F_00 = fidelity(ρ, P_00)

        push!(iterations, i)
        push!(fidelities, F_00)
        push!(successProbs, succProb)

    end

    distillable = false
    if F_00 > targetFid && i < maxIts
        distillable = true
    elseif isPpt || i >= maxIts
        distillable = false
    end

    return (distillable, iterations, fidelities, successProbs)

end

"""
    iterative_DEJMPS_protocol(ρ, targetFid, n, d, stdBasis, maxIts)

Applies iterations of the DEJMPS routine to the `n`-copy input state `ρ` in `d` dimensions. 
Iterates until `targetFid` or maximal number of iterations `maxIts` is reached. 
Returns `distillable=true` if `targetFid` could be reached and fidelities and success probabilies with respect to the `stdbasis` for each iteration.
"""
function iterative_DEJMPS_protocol(ρ, targetFid, n, d, stdbasis, maxIts=100)
    P_00 = stdbasis.basis[1][3]
    isPpt = isppt(ρ, d, 10)
    F_00 = fidelity(ρ, P_00)

    i = 0
    iterations = [i]
    fidelities = [F_00]
    successProbs = [0.0]

    while (F_00 <= targetFid) && !isPpt && (i < maxIts)

        i += 1

        (ρ, succProb) = DEJMPS_routine(ρ, n, d, stdbasis)
        isPpt = isppt(ρ, d, 10)

        F_00 = fidelity(ρ, P_00)

        push!(iterations, i)
        push!(fidelities, F_00)
        push!(successProbs, succProb)

    end

    distillable = false
    if F_00 > targetFid && i < maxIts
        distillable = true
    elseif isPpt || i >= maxIts
        distillable = false
    end

    return (distillable, iterations, fidelities, successProbs)

end

"""
    efficiency(distillable, iterations, successProbs, n)

Returns efficiency of `n`-copy distillation protocol with `iterations` iterations of success probabilities `successProbs`.
"""
function efficiency(distillable, iterations, successProbs, n)

    if !distillable
        return (0.0)
    else
        nIts = length(iterations[2:end])
        totalProb = prod(successProbs[2:end])

        if length(nIts) == 0
            eff = 0.0
        else
            eff = 1 / n^nIts * totalProb
        end
        return (eff)
    end

end