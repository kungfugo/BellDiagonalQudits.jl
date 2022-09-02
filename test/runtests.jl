using BellDiagonalQudits
using Test
using LazySets: tovrep, vertices_list
using LinearAlgebra: Diagonal

testStandardBasis3 = createStandardIndexBasis(3, 10)
testStandardKernel3 = createKernelPolytope(3, testStandardBasis3)
testExtKernel3 = extendSepVPolytopeBySepStates(tovrep(testStandardKernel3), [createDensityState(CoordState(1 / 9 * ones(9), "SEP"), testStandardBasis3)], 10)
testBipartiteWeylBasis3 = createBipartiteWeylOpereatorBasis(3)
testDictionaries3 = createDictionaryFromBasis(testStandardBasis3)
testMub3 = createStandardMub(3)
testRandomCoordWits = map(x -> getBoundedCoordEw(x), createRandomBoundedWits(3, testStandardBasis3, 2, true, 20))
testSyms3 = generateAllSymmetries(testStandardBasis3, 3)
testAnaSpec = AnalysisSpecification(true, true, true, true, true, true, true, false)
testAnaSpecSym = AnalysisSpecification(true, true, true, true, true, true, true, true)

@testset "BellDiagonalQudits.jl" begin

    # BellDiagonalQudits/src/BellStates
    @test length(createRandomCoordStates(10, 3, "magicSimplex")) == 10
    @test length(createStandardIndexBasis(4, 10).basis) == 16
    @test testDictionaries3[1][2, 1] == 6
    @test createDensityState(CoordState(1 / 9 * ones(9), "UNKNOWN"), testStandardBasis3).coords ≈ 1 / 9 * ones(9)
    @test length(createBipartiteWeylOpereatorBasis(3)) == 81

    # BellDiagonalQudits/src/EntanglementWitness
    @test length(testRandomCoordWits) == 2

    # BellDiagonalQudits/src/SeparableStates
    @test length(vertices_list(tovrep(testStandardKernel3))) == 12
    @test length(vertices_list(testExtKernel3)) == 13

    # BellDiagonalQudits/src/Symmetries
    @test length(testSyms3) == 432
    @test getSymCoords(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        generateAllSymmetries(testStandardBasis3, 3)[1:2]
    ) == [
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    ]

    # BellDiagonalQudits/src/EntanglementChecks
    @test all(
        map(
            i -> getfield(
                analyseCoordState(
                    3,
                    CoordState(1 / 9 * ones(9), "UNKNOWN"),
                    testAnaSpec,
                    testStandardBasis3,
                    testStandardKernel3,
                    testBipartiteWeylBasis3,
                    testDictionaries3,
                    testMub3,
                    testRandomCoordWits
                ), i
            ) == getfield(
                AnalysedCoordState(
                    CoordState(1 / 9 * ones(9), "UNKNOWN"),
                    true,
                    true,
                    true,
                    false,
                    false,
                    false,
                    false
                ), i),
            propertynames(
                AnalysedCoordState(
                    CoordState(1 / 9 * ones(9), "UNKNOWN"),
                    true,
                    true,
                    true,
                    false,
                    false,
                    false,
                    false
                )
            )[2:8]
        )
    )
    @test all(
        map(
            i -> getfield(
                symAnalyseCoordState(
                    3,
                    CoordState(1 / 9 * ones(9), "UNKNOWN"),
                    testSyms3[1:10],
                    testAnaSpecSym,
                    testStandardBasis3,
                    testStandardKernel3,
                    testBipartiteWeylBasis3,
                    testDictionaries3,
                    testMub3,
                    testRandomCoordWits
                ), i
            ) == getfield(
                AnalysedCoordState(
                    CoordState(1 / 9 * ones(9), "UNKNOWN"),
                    true,
                    true,
                    true,
                    false,
                    false,
                    false,
                    false
                ), i),
            propertynames(
                AnalysedCoordState(
                    CoordState(1 / 9 * ones(9), "UNKNOWN"),
                    true,
                    true,
                    true,
                    false,
                    false,
                    false,
                    false
                )
            )[2:8]
        )
    )
    @test map(
        x -> x.coordState.eClass,
        classifyAnalyzedStates!([
            AnalysedCoordState(
                CoordState(ones(9), "UNKNOWN"),
                true,
                true,
                true,
                false,
                false,
                false,
                false
            ),
            AnalysedCoordState(
                CoordState(ones(9), "UNKNOWN"),
                false,
                false,
                true,
                false,
                false,
                false,
                false
            ),
            AnalysedCoordState(
                CoordState(ones(9), "UNKNOWN"),
                false,
                false,
                true,
                true,
                false,
                false,
                false
            ),
            AnalysedCoordState(
                CoordState(ones(9), "UNKNOWN"),
                false,
                false,
                false,
                true,
                false,
                false,
                false
            )])
    ) == ["SEP", "PPT_UNKNOWN", "BOUND", "NPT"]
end
