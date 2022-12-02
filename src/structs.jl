"""
Represents a Bell basis related to Weyl operators.

- basis: Array with elements containing Bell basis density matrices, Weyl- and flat indices.
"""
struct StandardBasis
    basis::Array{Tuple{Int64,Tuple{Int64,Int64},Array{Complex{Float64},2}},1}
end;

"""
Represents a Bell diagonal state in Bell basis.

- coords: Coordinates in Bell basis. 
- eClass: Entanglement class of the represented state.
"""
mutable struct CoordState
    coords::Array{Float64,1}
    eClass::String
end;

"""
Represents a Bell diagonal state.

- coords: Coordinates in Bell basis.
- densityMatrix: Hermitian density matrix in computational basis.
- eClass: Entanglement class of the represented state.
"""
mutable struct DensityState
    coords::Array{Float64,1}
    densityMatrix::Hermitian{Complex{Float64},Array{Complex{Float64},2}}
    eClass::String
end;

"""
Represents an operator to detect Bell diagonal entangled states.

- coords: Coordinates in Bell basis
- operatorMatrix: Hermitian matrix representing the linear operator in computational basis.
"""
mutable struct EntanglementWitness
    coords::Array{Float64,1}
    operatorMatrix::Hermitian{Complex{Float64},Array{Complex{Float64},2}}
end;

"""
Represents an entanglement witness ``W`` with extrema and extremizers to detect entangled Bell diagonal states.

- coords: Coordinates in Bell basis.
- upperBound: Upper bound of ``tr W \\rho`` satisfied by all separable states ``\\rho``. Violation detects entanglement.
- lowerBound: Lower bound of ``tr W \\rho`` satisfied by all separable states ``\\rho``. Violation detects entanglement.
- maximizingDensityMatrix: Density matrix of separable state ``\\rho`` in computational basis, maximizing ``tr W \\rho``.
- minimizingDensityMatrix: Density matrix of separable state ``\\rho`` in computational basis. minimizing ``tr W \\rho``.
- checkedIterations: Number of iterations used in the optimization of bounds.
"""
mutable struct BoundedEW
    coords::Array{Float64,1}
    upperBound::Float64
    lowerBound::Float64
    maximizingDensityMatrix::Hermitian{Complex{Float64},Array{Complex{Float64},2}}
    minimizingDensityMatrix::Hermitian{Complex{Float64},Array{Complex{Float64},2}}
    checkedIterations::Int64
end;

"""
Represents an entanglement witness ``W`` with extrema to detect entangled Bell diagonal states.

- coords: Coordinates in Bell basis.
- upperBound: Upper bound of ``tr W \\rho`` satisfied by all separable states ``\\rho``. Violation detects entanglement.
- lowerBound: Lower bound of ``tr W \\rho`` satisfied by all separable states ``\\rho``. Violation detects entanglement.
- checkedIterations: Number of iterations used in the optimization of bounds.
"""
mutable struct BoundedCoordEW
    coords::Array{Float64,1}
    upperBound::Float64
    lowerBound::Float64
    checkedIterations::Int64
end;

"""
Specification which entanglement checks to use.

- kernel_check
- spinrepCheck
- ppt_check
- realignment_check 
- mub_check 
- numeric_ew_check
- useSymmetries
"""
mutable struct AnalysisSpecification
    kernel_check::Bool
    spinrepCheck::Bool
    ppt_check::Bool
    realignment_check::Bool
    concurrence_qp_check::Bool
    mub_check::Bool
    numeric_ew_check::Bool
    useSymmetries::Bool
end

"""
Representing an entanglement analyzed Bell diagonal state.

- coordState: The analyzed `CoordState`
- kernel: `true` if kernel check successful, else `false`. `missing` if entanglement check not applied.
- spinrep: `true` if spinrep check successful, else `false`. `missing` if entanglement check not applied.
- ppt: `true` if ppt check successful, else `false`. `missing` if entanglement check not applied.
- realign: `true` if realignment check successful, else `false`. `missing` if entanglement check not applied.
- concurrence: `true` if concurrence check successful, else `false`. `missing` if entanglement check not applied.
- mub: `true` if mub check successful, else `false`. `missing` if entanglement check not applied.
- numericEW: `true` if numericEW check successful, else `false`. `missing` if entanglement check not applied.
"""
mutable struct AnalysedCoordState
    coordState::CoordState
    kernel::Union{Missing,Bool}
    spinrep::Union{Missing,Bool}
    ppt::Union{Missing,Bool}
    realign::Union{Missing,Bool}
    concurrence::Union{Missing,Bool}
    mub::Union{Missing,Bool}
    numericEW::Union{Missing,Bool}
end;

"""
Exception for conflicts in analysis results.

- state: The `AnalysedCoordState` for which a conflict occurs.

"""
struct ClassConflictException <: Exception
    state::AnalysedCoordState
end;