struct StandardBasis
    basis::Array{Tuple{Int64,Tuple{Int64,Int64},Array{Complex{Float64},2}},1}
end;

mutable struct CoordState
    coords::Array{Float64,1}
    eClass::String
end;

mutable struct DensityState
    coords::Array{Float64,1}
    densityMatrix::Hermitian{Complex{Float64},Array{Complex{Float64},2}}
    eClass::String
end;

mutable struct EntanglementWitness
    coords::Array{Float64,1}
    operatorMatrix::Hermitian{Complex{Float64},Array{Complex{Float64},2}}
end;

mutable struct BoundedEW
    coords::Array{Float64,1}
    upperBound::Float64
    lowerBound::Float64
    maximizingDensityMatrix::Hermitian{Complex{Float64},Array{Complex{Float64},2}}
    minimizingDensityMatrix::Hermitian{Complex{Float64},Array{Complex{Float64},2}}
    checkedIterations::Int64
end;

mutable struct BoundedCoordEW
    coords::Array{Float64,1}
    upperBound::Float64
    lowerBound::Float64
    checkedIterations::Int64
end;

mutable struct AnalysisSpecification
    kernelCheck::Bool
    spinrepCheck::Bool
    pptCheck::Bool
    realignmentCheck::Bool
    concurrenceQpCheck::Bool
    mubCheck::Bool
    numericEwCheck::Bool
    useSymmetries::Bool
end

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

struct eClassConflictException <: Exception
    state::AnalysedCoordState
end;