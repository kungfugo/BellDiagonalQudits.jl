"""
    createKernelVertexStates(d, standardBasis::StandardBasis)

Return array containing collections of corresponding `standardBasis` indices, coordinates and density matrices for all `d` element sublattices in discrete phase space.
"""
function createKernelVertexStates(d, standardBasis::StandardBasis)

    allSubLattices = createDimElementSubLattices(d)
    vertexStates = map(x -> createIndexSubLatticeState(standardBasis, x), allSubLattices)

    return vertexStates
end

"""
    createKernelHPolytope(vertexCoordinates)

Return LazySets.HPolytope representation of polytope defined by `vertexCoordinates`.
"""
function createKernelHPolytope(vertexCoordinates)

    #Create polytope in vertex representation
    vPt = VPolytope(vertexCoordinates)

    #Convert to H-representation
    return tohrep(vPt)

end


"""
    createKernelPolytope(d, standardBasis::StandardBasis)

Return LazySets.HPolytope representation of the kernel polytope for dimension `d` and Bell basis `standardBasis`.
"""
function createKernelPolytope(d, standardBasis::StandardBasis)

    # Create kernel vertex states 
    vertexStates = createKernelVertexStates(d, standardBasis)

    # Get coordinates 
    vertexCoords = map(x -> x[2], vertexStates)

    # Create kernel polytope 
    return createKernelHPolytope(vertexCoords)

end

"""
    extendSepVPolytopeBySepStates(
        sepPolytope::VPolytope{Float64,Array{Float64,1}},
        sepDensityStates::Array{DensityState},
        precision::Integer
    )

Return an extended Lazysets.VPolytope representation of polytope of separable states based on given polytope `sepPolytope` and new separable `sepDensityStates` as new vertices.
"""
function extendSepVPolytopeBySepStates(
    sepPolytope::VPolytope{Float64,Array{Float64,1}},
    sepDensityStates::Array{DensityState},
    precision::Integer
)::VPolytope{Float64,Array{Float64,1}}

    existingVertices = vertices_list(sepPolytope)

    newVertices = unique(map(x -> roundDigits(abs.(x.coords), precision), sepDensityStates))

    combinedVertices = [existingVertices; newVertices]

    uniqueVertices = unique(map(x -> roundDigits(abs.(x), precision), combinedVertices))

    return VPolytope(uniqueVertices)

end