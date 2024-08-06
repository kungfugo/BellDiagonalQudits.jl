"""
    ispsd(M, precision)

Return `true` if the smallest eingenvalue of matrix `M` rounded to precision `precision` is not negative.
"""
function ispsd(M, precision)

    isPSD = false

    realEigvals = real.(eigvals(M))
    if round(minimum(realEigvals), digits=precision) >= 0
        isPSD = true
    end

    return isPSD

end

"""
    generic_vectorproduct(A,B)

For any vectors of equal length, return ``\\sum_i A[i]B[i]``, the sum of all products of elements with the same index.
"""
function generic_vectorproduct(A, B)

    @assert length(A) == length(B)

    d = length(A)

    Res = 0 * (A[1] * B[1])

    for i in (1:d)
        Res = Res + A[i] * B[i]
    end

    return Res
end

"""
    rounddigits(A, precision)

Return `A` with all elemets rounded up to `precision` digits.
"""
function rounddigits(A, precision)
    (x -> round(x, digits=precision)).(A)
end

"""
    isPPTP(ρ, d, precision)

Return `true` if the partial transposition of the Hermitian matrix based on the upper triangular of `ρ` is positive semi-definite in `precision`.
"""
function isppt(ρ, d, precision)

    M = Hermitian(ρ)
    minEigvalPT = real(ppt(M, [d, d], 2))

    return round(minEigvalPT, digits=precision) >= 0

end

"""
    create_dictionary_from_basis(stdBasis)

Return vector containing a dictionary and it's inverse, relating the ``d^2`` flat indices to the double indices of `stdBasis`.
"""
function create_dictionary_from_basis(stdBasis::StandardBasis)

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

"""
    map_indices_to_normcoords(indices, D)

Return `D` element normalized coordinate vector with equal nonzero values at given `indices`.
"""
function map_indices_to_normcoords(indices, D)

    c = zeros(D)
    for i in indices
        c[i] = 1 / length(indices)
    end

    return c

end

"""
    get_properdivisors(k:Int)

Return vector of proper divisors of `k`.
"""
function get_properdivisors(k::Int)

    divisors = Array{Int,1}(undef, 0)

    for j in (2:floor(sqrt(k)))
        if k % j == 0
            push!(divisors, j)
            if j != sqrt(k)
                push!(divisors, k / j)
            end
        end
    end

    return divisors
end


"""
    FT_OP(d)

Returns the `d`-dimensional Fourier gate.
"""
function FT_OP(d)
    F = complex(zeros(d, d))
    for i in 0:(d-1), j in (0:d-1)
        F += 1 / sqrt(d) * exp(2 * pi * im / d * j * i) * ketbra(i + 1, j + 1, d)
    end
    return (F)
end

"""
    GXOR_OP_ADD(d)

Returns the `d`-dimensional generalization of the GXOR gate acting as |i>k|j> => |i>|i+j>.
"""
function GXOR_OP_ADD(d)
    GXOR = complex(zeros(d^2, d^2))
    for i in 0:(d-1), m in 0:(d-1), j in 0:(d-1), n in 0:(d-1)
        GXOR += δ(i, j) * δ(m, mod(j + n, d)) * (modket(i, d) ⊗ modket(m, d)) * (modbra(j, d) ⊗ modbra(n, d))
    end
    return (GXOR)
end

"""
    GXOR_OP_SUB(d)

Returns the `d`-dimensional generalization of the GXOR gate acting as |i>k|j> => |i>|i-j>.
"""
function GXOR_OP_SUB(d)
    GXOR = complex(zeros(d^2, d^2))
    for i in 0:(d-1), m in 0:(d-1), j in 0:(d-1), n in 0:(d-1)
        GXOR += δ(i, j) * δ(m, mod(j - n, d)) * (modket(i, d) ⊗ modket(m, d)) * (modbra(j, d) ⊗ modbra(n, d))
    end
    return (GXOR)
end

"""
    depolarize_coords(coords)

Return depolarized elements of a coord probability vector. Leaves the first element invariant. Remaining elements are replaced by their average value.
"""
function depolarize_coords(coords)
    l = length(coords)
    newCoords = [coords[1]; (sum(coords[2:end]) / (l - 1)) * ones(l - 1)]
    return (newCoords)
end

"""
    belldiagonal_projection(stdbasis, ρ)

Projects any density matrix to the set of Bell-diagonal densitystates of corresponding dimension as defined by the `stdbasis`.
"""
function belldiagonal_projection(stdBasis::StandardBasis, ρ, eClass="UNKNOWN", precision::Integer=10)::DensityState
    coords = Float64[]
    for basisOp in stdBasis.basis
        push!(coords, real(round(tr(basisOp[3] * ρ), digits=precision)))
    end

    projState = create_densitystate(CoordState(coords, eClass), stdBasis)

    return projState
end

"""
    add_errorkeys(x,y,d,n)

Adds two error elements mod d.
"""
function add_errorkeys(x, y, d, n)

    @assert n == 2
    ((x_i1, x_j1), (x_i2, x_j2)) = (x[1], x[2])
    ((y_i1, y_j1), (y_i2, y_j2)) = (y[1], y[2])

    z_i1 = mod((x_i1 + y_i1), d)
    z_j1 = mod((x_j1 + y_j1), d)
    z_i2 = mod((x_i2 + y_i2), d)
    z_j2 = mod((x_j2 + y_j2), d)

    return (
        [(z_i1, z_j1), (z_i2, z_j2)]
    )

end

""" 
    modket(k,d)
Return the `k`th ket vector modulo `d`.
"""
function modket(k, d)
    return (ket(mod(k, d) + 1, d))
end

""" 
    modbra(k,d)
Return the `k`th bra vector modulo `d`.
"""
function modbra(k, d)
    return (bra(mod(k, d) + 1, d))
end

