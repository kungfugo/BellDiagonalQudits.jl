"""
    directOptimization(f, negf, method, d, iterations, constrainedOpt=false)

Return total (minimum, minimizer) and (maximum, maximizer) of `iterations` optimization runs of function `f` and its negative `negf` over the set of separale states using Optim.jl optimization mehtod `method`.

Optim.jl is used for optimazation based on `parameterization.jl`, so `f` and `negf` are defined for ``2(d-1)`` parameters. Supported methods include `NelderMead`, `LBFGS` and `NewtonTrustRegion`.
"""
function directOptimization(f, negf, method, d, iterations, constrainedOpt=false)

    lHalf = d - 1
    lMax = (2(d - 1))

    lowerBoundParam = zeros(2 * (d - 1))
    upperBoundParam = zeros(2 * (d - 1))

    for k in 1:lHalf
        upperBoundParam[k] = π / 2
    end
    for k in (lHalf+1):lMax
        upperBoundParam[k] = 2 * π
    end

    combinedLowerBounds = [lowerBoundParam; lowerBoundParam]
    combinedUpperBounds = [upperBoundParam; upperBoundParam]

    detMaxima = Float64[]
    detMinima = Float64[]
    detMaximizers = []
    detMinimizers = []

    for it in 1:iterations

        x_init_random = zeros(length(combinedUpperBounds))
        for i in eachindex(x_init_random)
            x_init_random[i] = rand(Uniform(combinedLowerBounds[i], combinedUpperBounds[i]))
        end

        if (!constrainedOpt)
            detOptmin = optimize(f, x_init_random, method())
            detOptmax = optimize(negf, x_init_random, method())
            detMin = Optim.minimum(detOptmin)
            detMax = -Optim.minimum(detOptmax)
        else
            detOptmin = optimize(f, combinedLowerBounds, combinedUpperBounds, x_init_random, Fminbox(method()))
            detOptmax = optimize(negf, combinedLowerBounds, combinedUpperBounds, x_init_random, Fminbox(method()))
            detMin = Optim.minimum(detOptmin)
            detMax = -Optim.minimum(detOptmax)
        end

        push!(detMaxima, detMax)
        push!(detMinima, detMin)
        push!(detMaximizers, detOptmax.minimizer)
        push!(detMinimizers, detOptmin.minimizer)

    end

    # Get total extrema of runs
    minimizerIndex = findmin(detMinima)
    maximizerIndex = findmax(detMaxima)
    totalMin = minimizerIndex[1]
    totalMax = maximizerIndex[1]
    totalMinimizer = detMinimizers[minimizerIndex[2]]
    totalMaximizer = detMaximizers[maximizerIndex[2]]

    return [
        (totalMin, totalMinimizer),
        (totalMax, totalMaximizer)
    ]

end