"""
    isPositiveSemiDefinite(M, precision)

Return `true`` if the smallest eingenvalue of matrix `M` rounded to precision `precision` is not negative.
"""
function isPositiveSemiDefinite(M, precision)

    isPSD = false

    realEigvals = real.(eigvals(M))
    if round(minimum(realEigvals), digits=precision) >= 0
        isPSD = true
    end

    return isPSD

end

"""
    genericVectorProduct(A,B)

For any vectors of equal length, return the sum of all products of elements with the same index.

Return ``\sum_i A[i]*B[i]``.
"""
function genericVectorProduct(A, B)

    @assert length(A) == length(B)

    d = length(A)

    Res = 0 * (A[1] * B[1])

    for i in (1:d)
        Res = Res + A[i] * B[i]
    end

    return Res
end

"""
    roundDigits(A, precision)

Return `A` with all elemets rounded up to `precision` digits.
"""
function roundDigits(A, precision)
    (x -> round(x, digits=precision)).(A)
end

"""
    isPPTP(ρ, d, precision)

Return `true` if the partial transposition of the Hermitian matrix based on the upper triangular of `ρ` is positive semi-definite within `precision` digits.
"""
function isPPT(ρ, d, precision)

    M = Hermitian(ρ)
    minEigvalPT = real(ppt(M, [d, d], 2))

    return round(minEigvalPT, digits=precision) >= 0

end

"""
    createDictionaryFromBasis(stdBasis)

Return vector containing a dictionarie and it's inverse, relating the d^2 flat indices to the double indices of a d-element `StandardBasis`.
"""
function createDictionaryFromBasis(stdBasis::StandardBasis)

    basisDict = Dict()
    for i in stdBasis.basis
        basisDict[i[2]] = i[1]
    end

    reverseDict = Dict()
    for i in stdBasis.basis
        reverseDict[i[1]] = i[2]
    end

    return [basisDict, reverseDict]
end

"""
    δ(n,m)

Return 1 if `n`==`m`
"""
function δ(n, m)::Int8
    if n == m
        return 1
    else
        return 0
    end
end

"""
    δ_mod(n,m,x)

Return 1 if `n` and `m` are congruent modulo `x`.
"""
function δ_mod(n, m, x)::Int8
    if mod(n, x) == mod(m, x)
        return 1
    else
        return 0
    end
end