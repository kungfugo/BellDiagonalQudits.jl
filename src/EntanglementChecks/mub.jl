# defining mubs

function createStandardMub(d)::Vector{Vector{Vector{ComplexF64}}}
    if d != 3
        throw("Only d=3 supported a this time.")
    end

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

    return (mubSet)

end

function calculateCorrelation(d, mubSet::Vector{Vector{Vector{ComplexF64}}}, ρ)

    C = 0
    for k in 1:length(mubSet)

        B = mubSet[k]
        if k == 1
            for i = 1:d
                C += tr(
                    Hermitian(proj(
                        B[i] ⊗ conj(B[mod((i - 1) + 2, d)+1])
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