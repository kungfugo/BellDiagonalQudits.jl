"""
kernel_check(coordState::CoordState, kernelPolytope::Union{HPolytope{Float64,Array{Float64,1}},VPolytope{Float64,Array{Float64,1}}})

Return `true` if the Euclidean coordinates of the `coordState`` are contained in the `kernelPolytope` represented in V- or H-representation.
"""
function kernel_check(coordState::CoordState, kernelPolytope::Union{HPolytope{Float64,Array{Float64,1}},VPolytope{Float64,Array{Float64,1}}})::Bool

    if coordState.coords ∈ kernelPolytope
        return true
    else
        return false
    end

end

"""
    ppt_check(coordState::CoordState, standardBasis::StandardBasis, precision=10)

Return `true`` if the `coordState` defined via the `standardBasis` has positive partial transposition in the given `precision`.
"""
function ppt_check(coordState::CoordState, standardBasis::StandardBasis, precision=10)::Bool

    densityState = create_densitystate(coordState, standardBasis)
    ρ = densityState.densityMatrix
    d = Int(sqrt(size(ρ, 1)))

    return isppt(ρ, d, precision)

end

"""
    realignment_check(coordState::CoordState, standardBasis::StandardBasis, precision=10)

Return `true`` if the realigned `coordState` defined via the `standardBasis` has trace norm ``> 1`` in the given `precision`.
"""
function realignment_check(coordState::CoordState, standardBasis::StandardBasis, precision=10)::Bool

    densityState = create_densitystate(coordState, standardBasis)
    ρ = densityState.densityMatrix

    r_ρ = reshuffle(ρ)

    return (round(norm_trace(r_ρ), digits=precision) > 1)

end

"""
    numeric_ew_check(coordState::CoordState, boundedEWs::Array{BoundedCoordEW}, relUncertainity::Float64)

Return `true` if any entanglement witness of `boundedEWs` detects the density matrix `ρ` as entangled.

An entanglement witness ``E`` of `boundedEWs` detects `ρ`, if the scalar product ``\\rho``.`coords` ``\\cdot E``.`coords` is not in [`lowerBound`, `upperBound`].
If a `relUncertainity` is given, the violation relative to `upperBound-lowerBound` needs to exceed `relUncertainity`` to detect entanglement. 
"""
function numeric_ew_check(coordState::CoordState, boundedEWs::Array{BoundedCoordEW}, relUncertainity=0.0)::Bool

    anyEntanglementFound = false

    for boundedEW in boundedEWs

        intervalLength = boundedEW.upperBound - boundedEW.lowerBound
        tolerance = intervalLength * relUncertainity

        anyEntanglementFound = !(
            (boundedEW.lowerBound - tolerance)
            <= dot(coordState.coords, boundedEW.coords)
            <= (boundedEW.upperBound + tolerance)
        )

        # If one EW witnesses entanglement, we can stop
        if anyEntanglementFound == true
            break
        end

    end

    return anyEntanglementFound

end

"""
    concurrence_qp_check(coordState::CoordState, d, dictionaries, precision=10)

Return `true` if the quasi-pure concurrence (see `concurrence.jl`) is positive for a `coordState` and given basis `dictionaries` in the given `precision`.
"""
function concurrence_qp_check(coordState::CoordState, d, dictionaries, precision=10)::Bool

    coords = coordState.coords
    if round(get_concurrence_qp(coords, d, dictionaries), digits=precision) > 0
        return true
    else
        return false
    end

end

"""
    mub_check(coordState::CoordState, d, stdBasis::StandardBasis, mubSet::Vector{Vector{Vector{ComplexF64}}})

Return `true` if the sum of mutual predictibilities for a `mubSet` (see `mub.jl`) of dimension `d` exceeds ``2`` for a `coordState` and given `standardBasis`.
"""
function mub_check(coordState::CoordState, d, stdBasis::StandardBasis, mubSet::Vector{Vector{Vector{ComplexF64}}})::Bool

    ρ = create_densitystate(coordState, stdBasis).densityMatrix

    if calculate_mub_correlation(d, mubSet, ρ) > 2
        return true
    else
        return false
    end

end

"""
    spinrep_check(coordState::CoordState, stdBasis::StandardBasis, bipartiteWeylBasis::Vector{Array{Complex{Float64},2}}, precision=10)

Return `true` and detects a `coordState` for a `standardBasis` as separbale, if its coefficiencts in the `bipartiteWeylBasis` have 1-norm smaller than ``2`` in given `precision`.
"""
function spinrep_check(coordState::CoordState, stdBasis::StandardBasis, bipartiteWeylBasis::Vector{Array{Complex{Float64},2}}, precision=10)

    ρ = create_densitystate(coordState, stdBasis).densityMatrix
    spinRepCoefficients = map(x -> tr(ρ * x'), bipartiteWeylBasis)

    absCoeffs = map(x -> real(sqrt(x' * x)), spinRepCoefficients)

    return (round(sum(absCoeffs), digits=precision) <= 2)

end
"""
    analyse_coordstate(
        d,
        coordState::CoordState,
        anaSpec::AnalysisSpecification,
        stdBasis::StandardBasis=missing,
        kernelPolytope::Union{HPolytope{Float64,Array{Float64,1}},VPolytope{Float64,Array{Float64,1}},Missing}=missing,
        bipartiteWeylBasis::Union{Vector{Array{Complex{Float64},2}},Missing}=missing,
        dictionaries::Union{Any,Missing}=missing,
        mubSet::Union{Vector{Vector{Vector{ComplexF64}}},Missing}=missing,
        boundedEWs::Union{Array{BoundedCoordEW},Missing}=missing,
        precision=10,
        relUncertainity=0.0
    )

Return an `AnalysedCoordState` for a `coordState` in `d` dimensions based on the given `anaSpec` and corresponding analysis objects.

If an entanglement check should not be carried out or if an analysis object in not passed as variable, the corresponding property in `anaSpec` needs to be `false`. 
In this case, return the corresponding property of the `AnalysedCoordState` as `missing`.
"""
function analyse_coordstate(
    d,
    coordState::CoordState,
    anaSpec::AnalysisSpecification,
    stdBasis::StandardBasis=missing,
    kernelPolytope::Union{HPolytope{Float64,Array{Float64,1}},VPolytope{Float64,Array{Float64,1}},Missing}=missing,
    bipartiteWeylBasis::Union{Vector{Array{Complex{Float64},2}},Missing}=missing,
    dictionaries::Union{Any,Missing}=missing,
    mubSet::Union{Vector{Vector{Vector{ComplexF64}}},Missing}=missing,
    boundedEWs::Union{Array{BoundedCoordEW},Missing}=missing,
    precision=10,
    relUncertainity=0.0
)::AnalysedCoordState

    anaCoordState = AnalysedCoordState(
        coordState,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing
    )

    # Kernel check
    if anaSpec.kernel_check && !ismissing(kernelPolytope)
        anaCoordState.kernel = kernel_check(coordState, kernelPolytope)
    end

    # Spinrep check
    if anaSpec.spinrepCheck && !ismissing(stdBasis) && !ismissing(bipartiteWeylBasis)
        anaCoordState.spinrep = spinrep_check(coordState, stdBasis, bipartiteWeylBasis, precision)
    end

    # PPT check
    if anaSpec.ppt_check && !ismissing(stdBasis)
        anaCoordState.ppt = ppt_check(coordState, stdBasis, precision)
    end

    # Realign check 
    if anaSpec.realignment_check && !ismissing(stdBasis)
        anaCoordState.realign = realignment_check(coordState, stdBasis, precision)
    end

    # Concurrence QP check 
    if anaSpec.concurrence_qp_check && !ismissing(dictionaries)
        anaCoordState.concurrence = concurrence_qp_check(coordState, d, dictionaries, precision)
    end

    # Mub check 
    if anaSpec.mub_check && !ismissing(stdBasis) && !ismissing(mubSet)
        anaCoordState.mub = mub_check(coordState, d, stdBasis, mubSet)
    end

    # numericEW check 
    if anaSpec.numeric_ew_check && !ismissing(boundedEWs)
        anaCoordState.numericEW = numeric_ew_check(coordState, boundedEWs, relUncertainity)
    end

    return anaCoordState

end

"""
    sym_analyse_coordstate(
        d,
        coordState::CoordState,
        symmetries::Array{Permutation},
        anaSpec::AnalysisSpecification,
        stdBasis::StandardBasis=missing,
        kernelPolytope::Union{HPolytope{Float64,Array{Float64,1}},VPolytope{Float64,Array{Float64,1}},Missing}=missing,
        bipartiteWeylBasis::Union{Vector{Array{Complex{Float64},2}},Missing}=missing,
        dictionaries::Union{Any,Missing}=missing,
        mubSet::Union{Vector{Vector{Vector{ComplexF64}}},Missing}=missing,
        boundedCoordEWs::Union{Array{BoundedCoordEW},Missing}=missing,
        precision=10,
        relUncertainity=0.0
    )

Return an `AnalysedCoordState` for a `coordState` in `d` dimensions based on the given `anaSpec` and corresponding analysis objects and symmetry analysis.

If an entanglement check should not be carried out or if an analysis object in not passed as variable, the corresponding property in `anaSpec` needs to be `false`. 
In this case, return the corresponding property of the `AnalysedCoordState` as `missing`.
"""
function sym_analyse_coordstate(
    d,
    coordState::CoordState,
    symmetries::Array{Permutation},
    anaSpec::AnalysisSpecification,
    stdBasis::StandardBasis=missing,
    kernelPolytope::Union{HPolytope{Float64,Array{Float64,1}},VPolytope{Float64,Array{Float64,1}},Missing}=missing,
    bipartiteWeylBasis::Union{Vector{Array{Complex{Float64},2}},Missing}=missing,
    dictionaries::Union{Any,Missing}=missing,
    mubSet::Union{Vector{Vector{Vector{ComplexF64}}},Missing}=missing,
    boundedCoordEWs::Union{Array{BoundedCoordEW},Missing}=missing,
    precision=10,
    relUncertainity=0.0
)::AnalysedCoordState

    if !anaSpec.useSymmetries
        throw("useSymmetries not specified in analysis specification")
    end

    copyAnaSpec = deepcopy(anaSpec)

    kernelCheckPassed = false
    spinrepCheckPassed = false
    pptCheckPassed = false
    realignmentCheckPassed = false
    concurrenceQpCheckPassed = false
    mubCheckPassed = false
    numericEwCheckPassed = false

    groupKernel = missing
    groupSpinrep = missing
    groupPpt = missing
    groupRealign = missing
    groupConcurrence = missing
    groupMub = missing
    groupNumericEw = missing

    # Create all symmetric states
    # Avoid duplicates 
    if length(coordState.coords) == length(unique(coordState.coords))
        symCoordStates = unique(map(
            x -> CoordState(x, coordState.eClass),
            get_symcoords(coordState.coords, symmetries)
        ))
    else
        symCoordStates = map(
            x -> CoordState(x, coordState.eClass),
            get_symcoords(coordState.coords, symmetries)
        )
    end

    # Analyse all symmetric states
    for symCoordState in symCoordStates

        analysedSymCoordState = analyse_coordstate(
            d,
            symCoordState,
            copyAnaSpec,
            stdBasis,
            kernelPolytope,
            bipartiteWeylBasis,
            dictionaries,
            mubSet,
            boundedCoordEWs,
            precision,
            relUncertainity
        )

        # Update copyAnaSpecs for sym group: Skip analysis for other group states if check was successful for state
        # Kernel check implies sep in kernel which is preserved ==> Keep searching if enabled and missing
        if copyAnaSpec.kernel_check
            kernelCheckPassed = !ismissing(analysedSymCoordState.kernel)
            if kernelCheckPassed
                groupKernel = analysedSymCoordState.kernel
            end
            copyAnaSpec.kernel_check = !kernelCheckPassed
        end

        # Spinrep check implies SEP ==> Keep searching if enabled and (false or missing)
        if copyAnaSpec.spinrepCheck
            spinrepCheckDone = !ismissing(analysedSymCoordState.spinrep)
            spinrepCheckPassed = !ismissing(analysedSymCoordState.spinrep) && analysedSymCoordState.spinrep
            if spinrepCheckPassed
                groupSpinrep = true
            elseif spinrepCheckDone
                groupSpinrep = false
            end
            copyAnaSpec.spinrepCheck = !spinrepCheckPassed
        end

        # Ppt check determines PPT/NPT which is preserved under symmetry ==> Keep searching if enabled and missing
        if copyAnaSpec.ppt_check
            pptCheckPassed = !ismissing(analysedSymCoordState.ppt)
            if pptCheckPassed
                groupPpt = analysedSymCoordState.ppt
            end
            copyAnaSpec.ppt_check = !pptCheckPassed
        end

        # Realignment check implies entanglement which is preserved under symmetry ==> Keep searching if enabled and (false or missing)
        if copyAnaSpec.realignment_check
            realignmentCheckDone = !ismissing(analysedSymCoordState.realign)
            realignmentCheckPassed = !ismissing(analysedSymCoordState.realign) && analysedSymCoordState.realign
            if realignmentCheckPassed
                groupRealign = true
            elseif realignmentCheckDone
                groupRealign = false
            end
            copyAnaSpec.realignment_check = !realignmentCheckPassed
        end

        # Concurrence check implies entanglement which is preserved under symmetry ==> Keep searching if enabled and (false or missing)
        if copyAnaSpec.concurrence_qp_check
            concurrenceQpCheckDone = !ismissing(analysedSymCoordState.concurrence)
            concurrenceQpCheckPassed = !ismissing(analysedSymCoordState.concurrence) && analysedSymCoordState.concurrence
            if concurrenceQpCheckPassed
                groupConcurrence = true
            elseif concurrenceQpCheckDone
                groupConcurrence = false
            end
            copyAnaSpec.concurrence_qp_check = !concurrenceQpCheckPassed
        end

        # Mub check implies entanglement which is preserved under symmetry ==> Keep searching if enabled and (false or missing)
        if copyAnaSpec.mub_check
            mubCheckDone = !ismissing(analysedSymCoordState.mub)
            mubCheckPassed = !ismissing(analysedSymCoordState.mub) && analysedSymCoordState.mub
            if mubCheckPassed
                groupMub = true
            elseif mubCheckDone
                groupMub = false
            end
            copyAnaSpec.mub_check = !mubCheckPassed
        end


        # EW check implies entanglement which is preserved under symmetry ==> Keep searching if enabled and (false or missing)
        if copyAnaSpec.numeric_ew_check
            numericEwCheckDone = !ismissing(analysedSymCoordState.numericEW)
            numericEwCheckPassed = !ismissing(analysedSymCoordState.numericEW) && analysedSymCoordState.numericEW
            if numericEwCheckPassed
                groupNumericEw = true
            elseif numericEwCheckDone
                groupNumericEw = false
            end
            copyAnaSpec.numeric_ew_check = !numericEwCheckPassed
        end

        allDetermined = !any([
            copyAnaSpec.kernel_check,
            copyAnaSpec.spinrepCheck,
            copyAnaSpec.ppt_check,
            copyAnaSpec.realignment_check,
            copyAnaSpec.concurrence_qp_check,
            copyAnaSpec.mub_check,
            copyAnaSpec.numeric_ew_check
        ])

        if allDetermined
            break
        end

    end

    anaSymCoordState = AnalysedCoordState(
        coordState,
        groupKernel,
        groupSpinrep,
        groupPpt,
        groupRealign,
        groupConcurrence,
        groupMub,
        groupNumericEw
    )

    return anaSymCoordState

end
"""
    classify_entanglement(analysedCoordState)

Return entanglement class of `analysedCoordState`. 

Entanglement class can be "UNKNWON", "PPT_UNKNOWN" for PPT states that can be separable or entangled, "SEP" for separable states, "BOUND" for PPT/bound entangled states or "NPT" for NPT/free entangled states.
"""
function classify_entanglement(analysedCoordState)

    class = "UNKNOWN"

    if (!ismissing(analysedCoordState.ppt)) && !analysedCoordState.ppt #npt
        class = "NPT"
    elseif (!ismissing(analysedCoordState.kernel) && analysedCoordState.kernel) || (!ismissing(analysedCoordState.spinrep) && analysedCoordState.spinrep) #separable
        class = "SEP"
    elseif (
        !ismissing(analysedCoordState.ppt)
        && analysedCoordState.ppt
        && !(
            (!ismissing(analysedCoordState.kernel) && analysedCoordState.kernel)
            ||
            (!ismissing(analysedCoordState.spinrep) && analysedCoordState.spinrep)
        )
    ) ##ppt not sep
        class = "PPT_UNKNOWN"

        if !ismissing(analysedCoordState.realign)
            if analysedCoordState.realign
                class = "BOUND"
            end
        end
        if !ismissing(analysedCoordState.concurrence)
            if analysedCoordState.concurrence
                class = "BOUND"
            end
        end
        if !ismissing(analysedCoordState.mub)
            if analysedCoordState.mub
                class = "BOUND"
            end
        end
        if !ismissing(analysedCoordState.numericEW)
            if analysedCoordState.numericEW
                class = "BOUND"
            end
        end
    end

    return class

end


"""
    classify_analyzed_states!(anaCoordStates::Array{AnalysedCoordState})

Set entanglement class for array of `analysedCoordStates`.
"""
function classify_analyzed_states!(analysedCoordStates::Array{AnalysedCoordState})
    for anaCoordState in analysedCoordStates

        derivedClass = classify_entanglement(anaCoordState)
        if derivedClass != "UNKNOWN"
            if anaCoordState.coordState.eClass == "UNKNOWN"
                anaCoordState.coordState.eClass = derivedClass
            else
                if anaCoordState.coordState.eClass != derivedClass
                    throw(ClassConflictException(anaCoordState))
                end
            end
        end
    end

    return analysedCoordStates
end