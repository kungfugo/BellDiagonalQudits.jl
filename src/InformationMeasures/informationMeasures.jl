"""
    coherent_information(ρ_AB, dims, subTrOut::Int)

Calculate the coherent information in bits.


# Arguments
- `ρ_AB`: Bipartite density matrix.
- `dims`: Tuple or vector of subsystem dimensions.
- `subTrOut`: The index of the subsystem to trace out (1 for subsystem A, 2 for B).

# Returns
- The coherent information value in bits (log base 2).
"""
function coherent_information(ρ_AB, dims, subTrOut::Int)
    
    ρ_AB = Hermitian(ρ_AB)
    subSystem = ptrace(Matrix(ρ_AB), dims, [subTrOut]) |> Hermitian

    jointEntropy = vonneumann_entropy(ρ_AB)
    subEntropy = vonneumann_entropy(subSystem)

    return((subEntropy - jointEntropy)/log(2))

end


"""
    mutual_information(ρ_AB, dims)

Calculate the Von Neumann mutual information in bits.

# Arguments
- `ρ_AB`: Bipartite density matrix.
- `dims`: Tuple or vector of subsystem dimensions `(d_A, d_B)`.

# Returns
- The mutual information value in bits (log base 2).
"""
function mutual_information(ρ_AB, dims)
    
    subSystem1 = ptrace(Matrix(ρ_AB), dims, [1]) |> Hermitian
    subSystem2 = ptrace(Matrix(ρ_AB), dims, [2]) |> Hermitian

    jointEntropy = vonneumann_entropy(ρ_AB)
    subEntropy1 = vonneumann_entropy(subSystem1)
    subEntropy2 = vonneumann_entropy(subSystem2)

    return((subEntropy1 + subEntropy2 - jointEntropy)/log(2))

end




""" 
    root_fidelity_sdp_primal(ρ, σ, returnVariables=false)

Calculate the root fidelity √F(ρ, σ) = tr|√ρ√σ| using a primal SDP formulation.

# Arguments
- `ρ`, `σ`: Density matrices of the same dimension.
- `returnVariables`: If `true`, also returns the optimal alignment matrix `X`.
"""
function root_fidelity_sdp_primal(ρ, σ, returnVariables=false)
    
    D = size(ρ, 1)

    model = Model(SCS.Optimizer)
    set_silent(model)

    @variable(model, X[1:D,1:D] in ComplexPlane())
    @objective(model, Max, real(tr(X)+tr(X')))
    @constraint(model, Hermitian([ρ X; X' σ]) in HermitianPSDCone())
    
    optimize!(model)

    if returnVariables
        return(1/2*objective_value(model),
            value.(X)
        )
    else
        return(1/2*objective_value(model) )
    end
    
end


 """ 
    root_fidelity_sdp_dual(ρ, σ, returnVariables=false)

Calculate the root fidelity √F(ρ, σ) using the dual SDP formulation.

# Arguments
- `ρ`, `σ`: Density matrices of the same dimension.
- `returnVariables`: If `true`, also returns the dual variables `Y` and `Z`.
"""   
function root_fidelity_sdp_dual(ρ, σ, returnVariables=false)
    
    D = size(ρ, 1)

    model = Model(SCS.Optimizer)
    set_silent(model)

    @variable(model, Y[1:D,1:D] in HermitianPSDCone())
    @variable(model, Z[1:D,1:D] in HermitianPSDCone())
    
    @objective(model, Min, real(tr(Y*ρ)+tr(Z*σ)))
    @constraint(model, Hermitian([Y I(D); I(D) Z]) in HermitianPSDCone())
    
    optimize!(model)

    if returnVariables
        return(1/2*objective_value(model),
            (value.(Y),value.(Z))
        )
    else
        return(1/2*objective_value(model) )
    end
    
end

""" 
    max_relative_entropy_sdp(ρ, σ, returnVariables=false)

Calculate the max-relative entropy.

# Arguments
- `ρ`, `σ`: Density matrices of the same dimension.
- `returnVariables`: If `true`, also returns the optimal scaling factor `λ`.
"""
function max_relative_entropy_sdp(ρ, σ, returnVariables=false)

    model = Model(SCS.Optimizer) 
    set_silent(model)
    
    @variable(model, λ)
    @objective(model, Min, λ)
    @constraint(model, λ >= 0)
    @constraint(model, Hermitian(λ*σ - ρ) in HermitianPSDCone())
    optimize!(model)

    if returnVariables
        return(objective_value(model) |> log2, value.(λ))
    else
        return(objective_value(model) |> log2)
    end
end

""" 
    max_relative_entropy_sdp_dual(ρ, σ, returnVariables=false)

Calculate the max-relative entropy using the dual SDP formulation.

# Arguments
- `ρ`, `σ`: Density matrices of the same dimension.
- `returnVariables`: If `true`, also returns the optimal witness operator `M`.
"""
function max_relative_entropy_sdp_dual(ρ, σ, returnVariables=false)

    d = size(ρ,1)
    @assert d == size(σ,1)

    model = Model(SCS.Optimizer) 
    set_silent(model)
    
    @variable(model, M[1:d, 1:d] in HermitianPSDCone())
    @objective(model, Max, real(tr(M*ρ)))
    @constraint(model, real(tr(M*σ))<=1)
    optimize!(model)
    

    if returnVariables
        return(objective_value(model) |> log2, value.(M))
    else
        return(objective_value(model) |> log2)
    end

end

""" 
    smooth_max_relative_entropy_primal(ρ, σ, ε, returnVariables=false)

Calculate the smooth max-relative entropy using the primal SDP.

# Arguments
- `ρ`, `σ`: Density matrices of the same dimension.
- `ε`: Smoothing parameter (fidelity-based).
- `returnVariables`: If `true`, returns the optimal SDP variables `(W, v, Z, u)`.
"""
function smooth_max_relative_entropy_primal(ρ, σ, ε, returnVariables=false)
    
    D = size(ρ, 1)

    model = Model(SCS.Optimizer)
    set_silent(model)

    @variable(model, W[1:D, 1:D] in HermitianPSDCone())
    @variable(model, v)
    @variable(model, Z[1:D, 1:D] in HermitianPSDCone())
    @variable(model, u)
    @objective(model, Max, u + 2*v*sqrt(1-ε) - real(tr(Z*ρ)))
    @constraint(model, v >= 0)
    @constraint(model, real(tr(W*σ)) <= 1)
    @constraint(model, Hermitian([Z v*I(D); v*I(D) (W-u*I(D))]) in HermitianPSDCone())
    
    optimize!(model)
    
    
    if returnVariables
        return(objective_value(model) |> log2, 
            (value.(W), value.(v), value.(Z), value.(u))
        )
    else
        return(objective_value(model) |> log2)
    end

end 

""" 
    smooth_max_relative_entropy_dual(ρ, σ, ε, returnVariables=false)

Calculate the smooth max-relative entropy using the dual SDP.

# Arguments
- `ρ`, `σ`: Density matrices of the same dimension.
- `ε`: Smoothing parameter (fidelity-based).
- `returnVariables`: If `true`, returns the optimal smoothed state `rho_AB` and alignment matrix `X`.
"""
function smooth_max_relative_entropy_dual(ρ, σ, ε, returnVariables=false)
    
    D = size(ρ, 1)

    model = Model(SCS.Optimizer)
    set_silent(model)

    
    @variable(model, λ >= 0)
    @variable(model, rho_AB[1:D, 1:D] in HermitianPSDCone())
    @variable(model, X[1:D,1:D] in ComplexPlane())
    
    @objective(model, Min, λ)
    @constraint(model, real(tr(rho_AB)) == 1)
    @constraint(model, Hermitian(λ*σ - rho_AB) in HermitianPSDCone())
    @constraint(model, real(tr(X)) >= sqrt(1-ε))
    @constraint(model, Hermitian([ρ X; X' rho_AB]) in HermitianPSDCone())
    optimize!(model)

    if returnVariables
        return(objective_value(model) |> log2,
            (value.(rho_AB),value.(X))
        )
    else
        return(objective_value(model) |> log2)
    end
    
end

""" 
    smooth_conditional_min_entropy_primal(ρ_AB, dims, ε, returnVariables=false)

Calculate the smooth conditional min-entropy using the primal SDP.

# Arguments
- `ρ_AB`: Bipartite density matrix.
- `dims`: Tuple `(d_A, d_B)` of subsystem dimensions.
- `ε`: Smoothing parameter.
- `returnVariables`: If `true`, returns the optimal SDP variables `(W, v, Z, u)`.
"""
function smooth_conditional_min_entropy_primal(ρ_AB, dims, ε, returnVariables=false)
    
    
    D = size(ρ_AB,1)
    (d_A, d_B) = dims
    ρ_A = Hermitian(partialtrace(ρ_AB, 2, [d_A,d_B]))
    ρ_B = Hermitian(partialtrace(ρ_AB, 1, [d_A,d_B]))

    model = Model(SCS.Optimizer)
    set_silent(model)

    @variable(model, W[1:D, 1:D] in HermitianPSDCone())
    @variable(model, v)
    @variable(model, Z[1:D, 1:D] in HermitianPSDCone())
    @variable(model, u)
    @objective(model, Max, -u + 2*v*sqrt(1-ε^2) - real(tr(Z*ρ_AB)))
    @constraint(model, v >= 0)
    @constraint(model, u >= 0)
    @constraint(model, Hermitian(I(d_B) - partialtrace(W, 1, [d_A, d_B])) in HermitianPSDCone())
    @constraint(model, Hermitian([Z v*I(D); v*I(D) (W + u*I(D))]) in HermitianPSDCone())
    
    optimize!(model)
    
    
    if returnVariables
        return(-(objective_value(model) |> log2), 
            (value.(W), value.(v), value.(Z), value.(u))
        )
    else
        return( -(objective_value(model) |> log2))
    end

end 

"""
    smooth_conditional_min_entropy_dual(ρ_AB, dims, ε, returnVariables=false)

Calculate the smooth conditional min-entropy using the dual SDP.


# Arguments
- `ρ_AB`: Bipartite density matrix.
- `dims`: Tuple `(d_A, d_B)` of subsystem dimensions.
- `ε`: Smoothing parameter (purified distance).
- `returnVariables`: If `true`, returns the optimal smoothed state and alignment matrix.
"""
function smooth_conditional_min_entropy_dual(ρ_AB, dims, ε, returnVariables=false)
    
    D = size(ρ_AB,1)
    (d_A, d_B) = dims
    ρ_A = Hermitian(partialtrace(ρ_AB, 2, [d_A,d_B]))
    ρ_B = Hermitian(partialtrace(ρ_AB, 1, [d_A,d_B]))

    model = Model(SCS.Optimizer)
    set_silent(model)

    
    @variable(model, rho_AB[1:D, 1:D] in HermitianPSDCone())
    @variable(model, S_B[1:d_B, 1:d_B] in HermitianPSDCone())
    @variable(model, X[1:D,1:D] in ComplexPlane())
    
    @objective(model, Min, real(tr(S_B)))
    @constraint(model, real(tr(rho_AB)) == 1) #Equal to <= Valid?
    @constraint(model, Hermitian(I(d_A)⊗S_B - rho_AB) in HermitianPSDCone())
    @constraint(model, real(tr(X)) >= sqrt(1-ε^2))
    @constraint(model, Hermitian([ρ_AB X; X' rho_AB]) in HermitianPSDCone())
    optimize!(model)

    if returnVariables
        return(-(objective_value(model) |> log2),
            (value.(rho_AB),value.(X))
        )
    else
        return(-(objective_value(model) |> log2))
    end
    
end

""" 
    smooth_conditional_max_entropy(ρ_AB, dims, ε, returnVariables=false)

Calculate the smooth conditional max-entropy via min-entropy duality.

# Arguments
- `ρ_AB`: Bipartite density matrix.
- `dims`: Tuple `(d_A, d_B)` of subsystem dimensions.
- `ε`: Smoothing parameter.
- `returnVariables`: If `true`, returns the variables from the underlying SDP.
"""
function smooth_conditional_max_entropy(ρ_AB, dims, ε, returnVariables=false)

    ψ_ABC = purify_density_matrix(Hermitian(ρ_AB))
    ρ_ABC = ψ_ABC*ψ_ABC'
    d_A = dims[1]
    d_B = dims[2]
    d_C = d_A*d_B

    ρ_AC = partialtrace(ρ_ABC, 2, [d_A, d_B, d_C])

    dims_AC = [d_A, d_C]
    
    smooth_cond_min_ent = smooth_conditional_min_entropy_primal(ρ_AC, dims_AC, ε, returnVariables)

    return(-smooth_cond_min_ent)

end


""" 
    max_mutual_information(ρ_AB, dims)

Calculate the Max-Mutual Information.

# Arguments
- `ρ_AB`: The bipartite density matrix.
- `dims`: Tuple or Vector `(d_A, d_B)` containing the dimensions of the subsystems.

Returns the mutual information in bits.
"""
function max_mutual_information(ρ_AB, dims)

    D = size(ρ_AB,1)
    (d_A, d_B) = dims
    ρ_A = Hermitian(partialtrace(ρ_AB, 2, [d_A,d_B]))
    ρ_B = Hermitian(partialtrace(ρ_AB, 1, [d_A,d_B]))

    return(max_relative_entropy_sdp(ρ_AB, ρ_A ⊗ ρ_B))

end

""" 
    smooth_max_mutual_information_routine_1_primal(ρ_AB, dims, ε, λ, precision=10, returnVariables=true, eps=1e-4, maxIters=100_000)

Execute "Routine 1" (Primal) from Popp et al. for smooth max-mutual information.

This SDP optimizes the smoothed state ρ' in the ε-ball of ρ_AB to maximize the 
feasibility parameter μ for a fixed λ.

# Arguments
- `ρ_AB`: Bipartite density matrix.
- `dims`: Dimensions (d_A, d_B).
- `ε`: Smoothing parameter (purified distance).
- `λ`: Fixed scaling factor for the marginal product ρ_A ⊗ ρ'_B.
- `eps`: Convergence tolerance for the SCS solver.
- `maxIters`: Maximum solver iterations.

Returns (μ, (ρ', X)) where μ > 0 implies λ is an upper bound for 2^I_max^ε.
"""
function smooth_max_mutual_information_routine_1_primal(ρ_AB, dims, ε, λ, precision=10, returnVariables=true, eps=1e-4, maxIters=100_000)

    D = size(ρ_AB,1)
    (d_A, d_B) = dims
    ρ_A = Hermitian(partialtrace(ρ_AB, 2, [d_A,d_B]))


    model = Model(SCS.Optimizer)
    set_optimizer_attribute(model, "eps_rel", eps) 
    set_optimizer_attribute(model, "max_iters", maxIters) 
    set_silent(model)
    @variable(model, μ >= 0)
    @variable(model, rho_AB[1:D, 1:D] in HermitianPSDCone())
    @variable(model, X[1:D,1:D] in ComplexPlane())
    @objective(model, Max, μ)
    @constraint(model, real(tr(rho_AB)) == 1)
    @constraint(model, Hermitian(λ*ρ_A⊗partialtrace(rho_AB, 1, [d_A,d_B]) - rho_AB - μ*I(D)) in HermitianPSDCone())
    @constraint(model, real(tr(X)) >= sqrt(1-ε^2))
    @constraint(model, Hermitian([ρ_AB X; X' rho_AB]) in HermitianPSDCone())
    optimize!(model)

    if returnVariables
        return(
            rounddigits(objective_value(model), precision),
            (value.(rho_AB),value.(X))
        )
    else
        return(rounddigits(objective_value(model), precision))
    end
    
end

""" 
    smooth_max_mutual_information_routine_1_dual(ρ_AB, dims, ε, λ, precision=10, returnVariables=true, eps=1e-4, maxIters=100_000)

Execute "Routine 1" (Dual) from Popp et al. to verify the feasibility of λ for I_max^ε.

This provides a dual certificate for the optimization. If the minimum objective value 
is non-negative, then λ is a valid upper bound for 2^I_max^ε.

# Arguments
- `ρ_AB`: Bipartite density matrix.
- `dims`: Tuple (d_A, d_B).
- `λ`: The candidate value for 2^I_max^ε.
"""
function smooth_max_mutual_information_routine_1_dual(ρ_AB, dims, ε, λ, precision=10, returnVariables=true, eps=1e-4, maxIters=100_000)

    D = size(ρ_AB,1)
    (d_A, d_B) = dims
    ρ_A = Hermitian(partialtrace(ρ_AB, 2, [d_A,d_B]))

    model = Model(SCS.Optimizer)
    set_optimizer_attribute(model, "eps_rel", eps) 
    set_optimizer_attribute(model, "max_iters", maxIters) 
    set_silent(model)
    @variable(model, n)
    @variable(model, v >= 0)
    @variable(model, G[1:D, 1:D] in HermitianPSDCone())
    @variable(model, W[1:D, 1:D] in HermitianPSDCone())
    @objective(model, Min, n - 2*v*sqrt(1-ε^2) + real(tr(G*ρ_AB)))
    @constraint(model, real(tr(W)) >= 1)
    @constraint(model, Hermitian([
        G v*I(D); 
        v*I(D) (n*I(D)-(λ*(I(d_A)⊗partialtrace(W*(ρ_A⊗I(d_B)),1,[d_A,d_B])) - W))
        ]) in HermitianPSDCone()
    )
    optimize!(model)

    if returnVariables
        return(
            rounddigits(objective_value(model), precision),
            (value.(n), value.(v), value.(G), value.(W))
        )
    else
        return(rounddigits(objective_value(model), precision))
    end
    
end


""" 
    smooth_max_mutual_information_routine_2(ρ_AB, dims, vrho_AB, precision=10)

Update the scaling factor λ for a fixed smoothed state `vrho_AB`.

This solves the SDP: min λ s.t. vrho_AB ≤ λ(ρ_A ⊗ vrho_B), where ρ_A is 
the marginal of the original state and vrho_B is the marginal of the smoothed state.

# Arguments
- `ρ_AB`: The original bipartite density matrix.
- `dims`: Tuple (d_A, d_B).
- `vrho_AB`: The smoothed density matrix found in Routine 1.
- `precision`: Rounding precision for the returned λ.
"""
function smooth_max_mutual_information_routine_2(ρ_AB, dims, vrho_AB, precision=10)

    D = size(ρ_AB,1)
    (d_A, d_B) = dims
    ρ_A = Hermitian(partialtrace(ρ_AB, 2, [d_A,d_B]))
    vrho_B = partialtrace(vrho_AB, 1, [d_A,d_B])

    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, λ >= 0)
    @objective(model, Min, λ)
    @constraint(model , Hermitian(λ*ρ_A⊗vrho_B - vrho_AB) in HermitianPSDCone())
    optimize!(model)

    return(rounddigits(objective_value(model), precision))
    
end

""" 
    smooth_max_mutual_information_approximation(ρ_AB, dims, ε, maxIts, approxOrder, precision=10, returnVariables=true)

Approximate the smooth max-mutual information using the iterative algorithm from Popp et al.

# Arguments
- `ρ_AB`: Bipartite density matrix.
- `dims`: Subsystem dimensions `(d_A, d_B)`.
- `ε`: Smoothing parameter.
- `maxIts`: Maximum number of see-saw iterations.
- `approxOrder`: Convergence threshold factor.
- `precision`: Rounding precision for SDP results.

Returns a tuple containing the final entropy (in bits), convergence status, and history of λ and μ.
"""
function smooth_max_mutual_information_approximation(ρ_AB, dims, ε, maxIts, approxOrder, precision=10, returnVariables=true)

    D = size(ρ_AB,1)
    (d_A, d_B) = dims

    ρ_A = Hermitian(partialtrace(ρ_AB, 2, [d_A,d_B]))
    ρ_B = Hermitian(partialtrace(ρ_AB, 1, [d_A,d_B]))
    
    λ_0 = 2^max_mutual_information(ρ_AB, dims)
    δ =  λ_0*float(10)^(-approxOrder)
    
    μ_results = []
    λ_results = [λ_0]

    λ = λ_0
    vrho_AB = ρ_AB

    converged = false

    for i in 1:maxIts

        (μ, (vrho_AB,vX)) =  smooth_max_mutual_information_routine_1_primal(ρ_AB, dims, ε, λ, precision, true)
        push!(μ_results, μ)
        if (μ <= 0 || isnan(μ))            
            @warn "Stopped loop due to small values of μ."
            break
        end

        λ = smooth_max_mutual_information_routine_2(ρ_AB, dims, vrho_AB, precision)
        push!(λ_results, λ)
        if (abs(λ_results[end]-λ_results[end-1]) <= δ)
            @info "Stopped loop due to convergence in approxOrder."
            converged = true
            break 
        end

        if i==maxIts 
            @warn "Reached maximal number of iterations."
        end

    end

    if returnVariables
        return(log2(λ), converged, λ_results, μ_results)
    else
        return(log2(λ), converged)
    end

end

""" 
    hypothesis_testing_relative_entropy(ρ, σ, ε)

Calculate the hypothesis testing relative entropy.

This implementation uses the dual SDP formulation for better numerical stability.

# Arguments
- `ρ`: The state corresponding to the null hypothesis.
- `σ`: The state corresponding to the alternative hypothesis.
- `ε`: The maximum allowed Type-I error (significance level).
"""
function hypothesis_testing_relative_entropy(ρ, σ, ε)

    d = size(ρ,1)
    @assert d == size(σ,1)
    
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, u)
    @variable(model, Z[1:d,1:d] in HermitianPSDCone())
    @objective(model, Max, real(u*(1-ε) - tr(Z)))
    @constraint(model, u >= 0)
    @constraint(model, Hermitian(σ + Z - u*ρ) in HermitianPSDCone())
   
    optimize!(model)

    return(-(objective_value(model) |> log2 ))
end

"""
    hypothesis_testing_mutual_information(ρ, dims, ε)

Calculate the hypothesis testing mutual information.

# Arguments
- `ρ`: Bipartite density matrix \$\\rho_{AB}\$.
- `dims`: A tuple or vector `(d_A, d_B)` specifying the dimensions of systems A and B.
- `ε`: The smoothing/error parameter (Type-I error tolerance), \$0 \\le \\varepsilon < 1\$.

# Returns
- The hypothesis testing mutual information in bits (log base 2).
"""
function hypothesis_testing_mutual_information(ρ, dims, ε)

    D = size(ρ,1)
    (d_A, d_B) = dims
    
    ρ_A = Hermitian(partialtrace(ρ, 2, [d_A,d_B]))
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, u)
    @variable(model, Z[1:D,1:D] in HermitianPSDCone())
    @variable(model, σ_B[1:d_B,1:d_B] in HermitianPSDCone())
    @objective(model, Max, real(u*(1-ε) - tr(Z)))
    @constraint(model, u >= 0)
    @constraint(model, real(tr(σ_B)) == 1)
    @constraint(model, Hermitian(ρ_A ⊗ σ_B + Z - u*ρ) in HermitianPSDCone())
   
    optimize!(model)

    return(-(objective_value(model) |> log2 ))
    
end

"""
    smooth_private_information(ρ, d, ε, ϵ, maxIts=50, approxOrder=5, precision=10)

Calculate the smooth_private_information of ρ.

The routine performs the following:
1. Purifies ρ to a state ψ_ABE.
2. Applies a local dephasing channel to system A, representing a measurement in the 
   computational basis.
3. Computes the hypothesis testing mutual information between Alice and Bob (ζ_AB).
4. Computes the smooth max-mutual information between Alice and Eve (ζ_AE).
5. Returns the difference: I_H^ε(A:B) - I_max^ϵ(A:E).

# Arguments
- `ρ`: Input density matrix (system A).
- `d`: Dimension of subsystem A and B (assumes d_A = d_B = d).
- `ε`: Error parameter for the hypothesis testing mutual information (Bob).
- `ϵ`: Smoothing parameter for the max-mutual information (Eve).
"""
function smooth_private_information(ρ, d, ε, ϵ, maxIts=50, approxOrder=5, precision=10)

    ψ_ABE = purify_density_matrix(Hermitian(ρ), precision)
    τ_ABE = ψ_ABE*ψ_ABE'

    ζ_ABE = local_dephasing_channel([d,d,d^2], 1)(τ_ABE)

    d_A, d_B, d_E = [d,d,d^2]

    ζ_AB = partialtrace(ζ_ABE, 3, [d_A, d_B, d_E])    
    ζ_AE = partialtrace(ζ_ABE, 2, [d_A, d_B, d_E])    

    hypoTestingMutualInf = hypothesis_testing_mutual_information(ζ_AB, [d_A, d_B], ε)
    smoothMaxMutualInf = smooth_max_mutual_information_approximation(ζ_AE, [d_A, d_E], ϵ, maxIts, approxOrder, precision, false)

    return hypoTestingMutualInf - smoothMaxMutualInf[1]

end
