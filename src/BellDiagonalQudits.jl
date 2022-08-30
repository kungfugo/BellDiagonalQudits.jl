module BellDiagonalQudits

using LinearAlgebra: eigvals, I, Hermitian, tr, dot, normalize, Diagonal
using QuantumInformation: proj, ket, ketbra, max_entangled, ⊗, reshuffle, norm_trace, ppt
using Distributions: Uniform, Normal
using Permutations: Permutation
using LazySets: VPolytope, HPolytope, tohrep, tovrep, vertices_list, ∈
using Polyhedra
using Optim

export
    CoordState, StandardBasis, DensityState, BoundedCoordEW, AnalysisSpecification, eClassConflictException, AnalysedCoordState,
    createRandomCoordStates, createStandardIndexBasis, createDensityState, createBipartiteWeylOpereatorBasis,
    createKernelPolytope, extendSepVPolytopeBySepStates,
    createRandomBoundedWits, getBoundedCoordEw,
    generateAllSymmetries, getSymCoords,
    analyseCoordState, symAnalyseCoordState, classifyAnalyzedStates!,
    createDictionaryFromBasis, createStandardMub

include("structs.jl")
include("Utils/utils.jl")
include("BellStates/bellStates.jl")
include("Parameterization/parameterization.jl")
include("Optimization/optimization.jl")
include("EntanglementWitnesses/eWitnesses.jl")
include("Symmetries/symmetries.jl")
include("EntanglementChecks/concurrence.jl")
include("EntanglementChecks/mub.jl")
include("EntanglementChecks/entanglementChecks.jl")
include("SeparableStates/separableStates.jl")

end


