function getPermutationFromTranslation(translation::Tuple{Int,Int}, stdBasis::StandardBasis, d)

    myDicts = createDictionaryFromBaxsis(stdBasis)
    dict = myDicts[1] #keys are tuples 
    reverseDict = myDicts[2] #keys are indices

    D = length(myDicts[1])
    t = collect(translation)

    perm = zeros(D)

    for i in 1:D

        oldCoords = collect(reverseDict[i])
        newCoords = mod.(oldCoords + t, d)
        permIndex = dict[Tuple(newCoords)]

        perm[i] = permIndex

    end

    return Permutation(perm)

end

function getPermutationFromMomentumInversion(stdBasis::StandardBasis, d)

    #load dictionary
    myDicts = createDictionaryFromBaxsis(stdBasis)
    dict = myDicts[1] #keys are tuples 
    reverseDict = myDicts[2] #keys are indices

    D = length(myDicts[1])

    perm = zeros(D)

    for i in 1:D

        oldCoords = collect(reverseDict[i])
        newCoords = mod.([-oldCoords[1], oldCoords[2]], d)
        permIndex = dict[Tuple(newCoords)]

        perm[i] = permIndex

    end

    return Permutation(perm)

end

function getPermutationFromQuarterRotation(stdBasis::StandardBasis, d)

    #load dictionary
    myDicts = createDictionaryFromBaxsis(stdBasis)
    dict = myDicts[1] #keys are tuples 
    reverseDict = myDicts[2] #keys are indices

    D = length(myDicts[1])

    perm = zeros(D)

    for i in 1:D

        oldCoords = collect(reverseDict[i])
        newCoords = mod.([oldCoords[2], -oldCoords[1]], d)
        permIndex = dict[Tuple(newCoords)]

        perm[i] = permIndex

    end

    return Permutation(perm)

end

function getPermutationFromVerticalShear(stdBasis::StandardBasis, d)

    #load dictionary
    myDicts = createDictionaryFromBaxsis(stdBasis)
    dict = myDicts[1] #keys are tuples 
    reverseDict = myDicts[2] #keys are indices

    D = length(myDicts[1])

    perm = zeros(D)

    for i in 1:D

        oldCoords = collect(reverseDict[i])
        newCoords = mod.([oldCoords[1] + oldCoords[2], oldCoords[2]], d)
        permIndex = dict[Tuple(newCoords)]

        perm[i] = permIndex

    end

    return Permutation(perm)

end

function generateAllSymmetries(stdBasis::StandardBasis, d, orderLimit=0)::Array{Permutation}

    allSymmetries = Permutation[]
    #Get all translations
    for i in 0:(d-1)
        for j in 0:(d-1)
            push!(allSymmetries, getPermutationFromTranslation((i, j), stdBasis, d))
        end
    end

    generators = (
        getPermutationFromQuarterRotation(stdBasis, d),
        getPermutationFromMomentumInversion(stdBasis, d),
        getPermutationFromVerticalShear(stdBasis, d)
    )

    findNew = true
    order = 1

    while findNew

        foundSims = Permutation[]

        for x in generators
            generatedSyms = map(s -> s * x, allSymmetries)
            newSyms = filter(e -> !(e in allSymmetries), generatedSyms)
            append!(foundSims, newSyms)
        end

        if length(foundSims) > 0
            append!(allSymmetries, foundSims)
        end

        order += 1

        findNew = (length(foundSims) > 0 && (orderLimit >= order || orderLimit == 0))

        unique!(allSymmetries)

    end

    return allSymmetries

end

function getSymCoords(coords::Array{Float64,1}, symPermutations::Array{Permutation})::Array{Array{Float64,1}}

    symCoords = map(x -> Matrix(x) * coords, symPermutations)
    return symCoords

end