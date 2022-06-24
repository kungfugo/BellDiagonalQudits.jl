function isPositiveSemiDefinite(M, precision)

    isPSD = false

    realEigvals = real.(eigvals(M))
    if round(minimum(realEigvals), digits=precision) >= 0
        isPSD = true
    end

    return isPSD

end

function myMulti(coord, iB)
    d = length(iB)
    Res = zeros(d, d)

    for i in (1:d)
        Res = Res + coord[i] * iB[i]
    end

    return Res
end

function roundDigits(A, b)
    (x -> round(x, digits=b)).(A)
end

function isPPT(ρ, d, precision)

    M = Hermitian(ρ)
    minEigvalPT = real(ppt(M, [d, d], 2))

    return round(minEigvalPT, digits=precision) >= 0

end

function createDictionaryFromBaxsis(stdBasis::StandardBasis)

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

#discrete delta function 
function δ(n, m)::Int8
    if n == m
        return 1
    else
        return 0
    end
end

function δ_mod(n, m, x)::Int8
    if mod(n, x) == mod(m, x)
        return 1
    else
        return 0
    end
end