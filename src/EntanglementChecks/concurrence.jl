"""
    get_concurrence_qp(coords, d, dictionaries)

Return quasi-pure approximation of the concurrence for a Bell diagonal state represented by coordinates `coords` with respect to a `StandardBasis` and corresponding `dictionaries` in `d` dimensions.
"""
function get_concurrence_qp(coords, d, dictionaries)::Float64

    S = Float64[]
    D = length(coords)
    dict = dictionaries[2]
    revDict = dictionaries[1]

    (maxCoord, maxIndex) = findmax(coords)
    (n, m) = dict[maxIndex]

    for i in 1:D

        (k, l) = dict[i]
        push!(
            S,
            sqrt(
                max(
                    0,
                    d / (2 * (d - 1)) * coords[i] * (
                        (1 - 2 / d) * maxCoord * δ(i, maxIndex)
                        +
                        1 / d^2 * coords[revDict[(
                            mod(2 * n - k, d),
                            mod(2 * m - l, d)
                        )]]
                    )
                )
            )
        )
    end

    maxS = findmax(S)[1]
    cQP = maxS - (sum(S) - maxS)
    if cQP > 0
        return cQP
    else
        return 0
    end
end

"""
    createProjectorOperator(d)

Return concurrence related operator in arbitrary dimension `d`.
"""
function createProjectorOperator(d)

    A = proj(0 * ket(1, d^4))

    for j in 1:d
        for i in 1:(j-1)
            for l in 1:d
                for k in 1:(l-1)

                    A += proj(
                        ket(i, d) ⊗ ket(k, d) ⊗ ket(j, d) ⊗ ket(l, d)
                        -
                        ket(j, d) ⊗ ket(k, d) ⊗ ket(i, d) ⊗ ket(l, d)
                        -
                        ket(i, d) ⊗ ket(l, d) ⊗ ket(j, d) ⊗ ket(k, d)
                        +
                        ket(j, d) ⊗ ket(l, d) ⊗ ket(i, d) ⊗ ket(k, d)
                    )
                end
            end
        end
    end

    return 4 * A

end

""" 
    get_concurrence_qp_gendiagonal(coords, d, basisStates)

Return quasi-pure approximation of the concurrence for a `d`-dimensional Bell diagonal state represented by coordinates `coords` with respect to a set of `basisStates`. 
"""
function get_concurrence_qp_gendiagonal(coords, d, basisStates::Array{Vector{ComplexF64}})

    maxIndex = findmax(coords)[2]
    A = createProjectorOperator(d)
    Ω_max = basisStates[maxIndex]
    Χ = A * (Ω_max ⊗ Ω_max)
    Χ = Χ / norm(Χ)


    T = complex(zeros(d^2, d^2))

    for i in 1:d^2, j in 1:d^2

        μ_i = coords[i]
        μ_j = coords[j]

        T[i, j] = sqrt(μ_i * μ_j) * (basisStates[i]' ⊗ basisStates[j]') * Χ

    end

    sV = svd(T).S
    maxSv = findmax(sV)[1]

    concurrence_qp = max(0, 2 * maxSv - sum(sV))

    return (concurrence_qp)

end