module MagicSimplex

using LinearAlgebra: eigvals, I, Hermitian, tr, dot
using QuantumInformation: proj, ket, ketbra, max_entangled, ⊗, reshuffle, norm_trace, ppt
using Distributions: Uniform
using Permutations: Permutation
using LazySets: HPolytope, VPolytope, ∈
using JLD: load, save

export
    CoordState, StandardBasis, DensityState, BoundedCoordEW, AnalysisSpecification, eClassConflictException, AnalysedCoordState,
    createRandomCoordStates, createStandardIndexBasis, createDensityState,
    generateAllSymmetries, getSymCoords,
    analyseCoordState, symAnalyseCoordState, classifyAnalyzedStates!,
    load, save,
    createDictionaryFromBaxsis, createStandardMub

include("structs.jl")
include("Utils/utils.jl")
include("BellStates/bellStates.jl")
include("Symmetries/symmetries.jl")
include("EntanglementChecks/concurrence.jl")
include("EntanglementChecks/mub.jl")
include("EntanglementChecks/entanglementChecks.jl")

end


