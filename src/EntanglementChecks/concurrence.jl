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