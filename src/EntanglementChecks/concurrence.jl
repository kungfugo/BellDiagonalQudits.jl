# Implementing analyses related to the 'concurrence' as entanglement releated function.
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

##Explicit formula
function getConcurrenceQP(coords, d, dictionaries)::Float64

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