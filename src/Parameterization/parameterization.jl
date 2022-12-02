"""
    create_red_parmatrix_from_parvector(x, d)

Return parameter matrix for pure state parameterization from parameter vector `x` of length ``2(d-1)``.
"""
function create_red_parmatrix_from_parvector(x, d)

    if length(x) != (2 * (d - 1))
        throw("Param vector and dimension mismatch")
    end

    λ = zeros(d, d)
    k = 0

    #fill first row with first elements of x
    for n in 2:d
        k = k + 1
        λ[1, n] = x[k]
    end

    # fill first column with second half elements of x
    for n in 2:d
        k = k + 1
        λ[n, 1] = x[k]
    end

    return λ

end


"""
    create_comppar_unitary(λ, d)

Return `d` dimensional prameterized unitary matrix from parameter matrix `λ`.
"""
function create_comppar_unitary(λ, d)

    # create onb
    P = Array{
        Array{Complex{Float64},2},
        1}(undef, 0)

    for i in (1:d)
        b_i = proj(ket(i, d))
        push!(P, b_i)
    end

    σ = Array{Array{Complex{Float64},2},2}(undef, d, d)


    # define factor 1 for global phases 
    R = I(d)
    for l in 1:d
        R = R * exp(im * P[l] * λ[l, l])
    end

    # define factor 2 for rotations and relative phases
    L = I(d)
    for m in (1:d-1)
        for n in (m+1:d)
            σ[m, n] = -im * ketbra(m, n, d, d) + im * ketbra(n, m, d, d)
            L = L * exp(im * P[n] * λ[n, m]) * exp(im * σ[m, n] * λ[m, n])
        end
    end

    U_C = L * R

    return U_C

end

"""
    get_comppar_unitary_from_parvector(x, d)

Return parameterized unitatry matrix U of dimension `d` and rank 1 from parameter vector `x` with ``2(d-1)`` elements. 

Using the first basis state of the computational basis with density matrix ``e_1``, any pure state ``\\rho`` can be generated as ``\\rho = U e_1 U^\\dagger``.
"""
function get_comppar_unitary_from_parvector(x, d)
    λ = create_red_parmatrix_from_parvector(x, d)
    U = create_comppar_unitary(λ, d)
    return U
end