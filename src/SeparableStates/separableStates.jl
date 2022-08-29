"""
    createKernelVertexStates(d, standardBasis::StandardBasis)

Return array containing collections of corresponding `standardBasis` indices, coordinates and density matrix for all `d` element sublattices in discrete phase space.
"""
function createKernelVertexStates(d, standardBasis::StandardBasis)

    allSubLattices = createDimElementSubLattices(d)
    vertexStates = map(x -> createIndexSubLatticeState(standardBasis, x), allSubLattices)

    return vertexStates
end

"""
    createKernelHPolytope(vertexCoordinates)

Return LazySets: HPolytope defined by `vertexCoordinates`.
"""
function createKernelHPolytope(vertexCoordinates)

    #Create polytope in vertex representation
    vPt = VPolytope(vertexCoordinates)

    #Convert to H-representation
    return tohrep(vPt)

end


"""
    createKernelPolytope(d, standardBasis::StandardBasis)

Return LazySets: HPolytope corresponding the kernel polytope for dimension `d`.
"""
function createKernelPolytope(d, standardBasis::StandardBasis)

    # Create kernel vertex states 
    vertexStates = createKernelVertexStates(d, standardBasis)

    # Get coordinates 
    vertexCoords = map(x -> x[2], vertexStates)

    # Create kernel polytope 
    return createKernelHPolytope(vertexCoords)

end