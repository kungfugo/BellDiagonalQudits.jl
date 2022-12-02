"""
    create_kernel_vertexstates(d, standardBasis::StandardBasis)

Return array containing collections of corresponding `standardBasis` indices, coordinates and density matrices for all `d` element sublattices in discrete phase space.
"""
function create_kernel_vertexstates(d, standardBasis::StandardBasis)

    allSubLattices = create_dimelement_sublattices(d)
    vertexStates = map(x -> create_index_sublattice_state(standardBasis, x), allSubLattices)

    return vertexStates
end

"""
    create_kernel_hpolytope(vertexCoordinates)

Return LazySets.HPolytope representation of polytope defined by `vertexCoordinates`.
"""
function create_kernel_hpolytope(vertexCoordinates)

    #Create polytope in vertex representation
    vPt = VPolytope(vertexCoordinates)

    #Convert to H-representation
    return tohrep(vPt)

end


"""
    create_kernel_polytope(d, standardBasis::StandardBasis)

Return LazySets.HPolytope representation of the kernel polytope for dimension `d` and Bell basis `standardBasis`.
"""
function create_kernel_polytope(d, standardBasis::StandardBasis)

    # Create kernel vertex states 
    vertexStates = create_kernel_vertexstates(d, standardBasis)

    # Get coordinates 
    vertexCoords = map(x -> x[2], vertexStates)

    # Create kernel polytope 
    return create_kernel_hpolytope(vertexCoords)

end

"""
    extend_vpolytope_by_densitystates(
        sepPolytope::VPolytope{Float64,Array{Float64,1}},
        sepDensityStates::Array{DensityState},
        precision::Integer
    )

Return an extended Lazysets.VPolytope representation of polytope of separable states based on given polytope `sepPolytope` and new separable `sepDensityStates` as new vertices.
"""
function extend_vpolytope_by_densitystates(
    sepPolytope::VPolytope{Float64,Array{Float64,1}},
    sepDensityStates::Array{DensityState},
    precision::Integer
)::VPolytope{Float64,Array{Float64,1}}

    existingVertices = vertices_list(sepPolytope)

    newVertices = unique(map(x -> rounddigits(abs.(x.coords), precision), sepDensityStates))

    combinedVertices = [existingVertices; newVertices]

    uniqueVertices = unique(map(x -> rounddigits(abs.(x), precision), combinedVertices))

    return VPolytope(uniqueVertices)

end