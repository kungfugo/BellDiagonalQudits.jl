module BellDiagonalQudits

using LinearAlgebra: eigvals, I, Hermitian, tr, dot, normalize, Diagonal
using QuantumInformation: proj, ket, ketbra, max_entangled, ⊗, reshuffle, norm_trace, ppt
using Distributions: Uniform, Normal
using Permutations: Permutation
using LazySets: VPolytope, HPolytope, tohrep, tovrep, vertices_list, ∈
using Polyhedra
using Optim

export
    CoordState, StandardBasis, DensityState, BoundedCoordEW, AnalysisSpecification, ClassConflictException, AnalysedCoordState,
    uniform_bell_sampler, create_random_coordstates, create_standard_indexbasis, create_densitystate, create_bipartite_weyloperator_basis,
    create_kernel_polytope, extend_vpolytope_by_densitystates,
    create_random_bounded_ews, get_bounded_coordew,
    generate_symmetries, get_symcoords,
    analyse_coordstate, sym_analyse_coordstate, classify_analyzed_states!,
    create_dictionary_from_basis, create_standard_mub

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


