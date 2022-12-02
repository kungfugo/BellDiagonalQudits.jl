"""
    create_standard_mub(d)

Return vector of mutually unbiased bases for dimensions `d` three or four.
"""
function create_standard_mub(d)::Vector{Vector{Vector{ComplexF64}}}

    mubSet = Vector{Vector{ComplexF64}}[]

    if d == 3
        #There are 4 mutually unbiased Bases. 
        mubSet = Vector{Vector{ComplexF64}}[]

        w = exp(2 // 3 * π * im)
        B1 = [
            ket(1, 3),
            ket(2, 3),
            ket(3, 3)
        ]
        B2 = [
            1 / sqrt(3) * (ket(1, 3) + ket(2, 3) + ket(3, 3)),
            1 / sqrt(3) * (ket(1, 3) + w * ket(2, 3) + w^2 * ket(3, 3)),
            1 / sqrt(3) * (ket(1, 3) + w^2 * ket(2, 3) + w * ket(3, 3))
        ]
        B3 = [
            1 / sqrt(3) * (ket(1, 3) + w * ket(2, 3) + w * ket(3, 3)),
            1 / sqrt(3) * (ket(1, 3) + w^2 * ket(2, 3) + ket(3, 3)),
            1 / sqrt(3) * (ket(1, 3) + ket(2, 3) + w^2 * ket(3, 3))
        ]
        B4 = [
            1 / sqrt(3) * (ket(1, 3) + w^2 * ket(2, 3) + w^2 * ket(3, 3)),
            1 / sqrt(3) * (ket(1, 3) + ket(2, 3) + w * ket(3, 3)),
            1 / sqrt(3) * (ket(1, 3) + w * ket(2, 3) + ket(3, 3))
        ]

        push!(mubSet, B1)
        push!(mubSet, B2)
        push!(mubSet, B3)
        push!(mubSet, B4)

    elseif d == 4

        B1 = [
            ket(1, 4),
            ket(2, 4),
            ket(3, 4),
            ket(4, 4)
        ]
        B2 = [
            1 / 2 * (ket(1, 4) + ket(2, 4) + ket(3, 4) + ket(4, 4)),
            1 / 2 * (ket(1, 4) + ket(2, 4) - ket(3, 4) - ket(4, 4)),
            1 / 2 * (ket(1, 4) - ket(2, 4) - ket(3, 4) + ket(4, 4)),
            1 / 2 * (ket(1, 4) - ket(2, 4) + ket(3, 4) - ket(4, 4))
        ]
        B3 = [
            1 / 2 * (ket(1, 4) - im * ket(2, 4) - im * ket(3, 4) - ket(4, 4)),
            1 / 2 * (ket(1, 4) - im * ket(2, 4) + im * ket(3, 4) + ket(4, 4)),
            1 / 2 * (ket(1, 4) + im * ket(2, 4) + im * ket(3, 4) - ket(4, 4)),
            1 / 2 * (ket(1, 4) + im * ket(2, 4) - im * ket(3, 4) + ket(4, 4))
        ]
        B4 = [
            1 / 2 * (ket(1, 4) - ket(2, 4) - im * ket(3, 4) - im * ket(4, 4)),
            1 / 2 * (ket(1, 4) - ket(2, 4) + im * ket(3, 4) + im * ket(4, 4)),
            1 / 2 * (ket(1, 4) + ket(2, 4) + im * ket(3, 4) - im * ket(4, 4)),
            1 / 2 * (ket(1, 4) + ket(2, 4) - im * ket(3, 4) + im * ket(4, 4))
        ]
        B5 = [
            1 / 2 * (ket(1, 4) - im * ket(2, 4) - ket(3, 4) - im * ket(4, 4)),
            1 / 2 * (ket(1, 4) - im * ket(2, 4) + ket(3, 4) + im * ket(4, 4)),
            1 / 2 * (ket(1, 4) + im * ket(2, 4) + ket(3, 4) - im * ket(4, 4)),
            1 / 2 * (ket(1, 4) + im * ket(2, 4) - ket(3, 4) + im * ket(4, 4))
        ]

        push!(mubSet, B1)
        push!(mubSet, B2)
        push!(mubSet, B3)
        push!(mubSet, B4)
        push!(mubSet, B5)

    else
        throw("Only d=3 and d=4 are supported")
    end

    return (mubSet)

end

"""
    calculate_mub_correlation(d, mubSet::Vector{Vector{Vector{ComplexF64}}}, ρ, s=-1)

Based on complete set of mutually unbiased bases `mubSet`, return sum of mutual predictibilities, shifted by `s`, for density matrix `ρ` in `d` dimensions. 
"""
function calculate_mub_correlation(d, mubSet::Vector{Vector{Vector{ComplexF64}}}, ρ, s=-1)

    if s == -1
        if d == 3
            s = 2
        elseif d == 4
            s = 3
        else
            s = 0
        end
    end

    C = 0
    for k in 1:length(mubSet)

        B = mubSet[k]
        if k == 1
            for i = 1:d
                C += tr(
                    Hermitian(proj(
                        B[i] ⊗ B[mod((i - 1) + s, d)+1]
                    ) * ρ)
                )
            end
        else
            for i = 1:d
                C += tr(
                    Hermitian(proj(
                        B[i] ⊗ conj(B[i])
                    ) * ρ)
                )
            end
        end
    end

    return (C)

end