using BellDiagonalQudits
using Test
using LazySets: tovrep, vertices_list
using JLD: load

testStandardBasis3 = createStandardIndexBasis(3, 10)
testStandardKernel3 = createKernelPolytope(3, testStandardBasis3)
testBipartiteWeylBasis3 = createBipartiteWeylOpereatorBasis(3)
testDictionaries3 = createDictionaryFromBasis(testStandardBasis3)
testMub3 = createStandardMub(3)
testAnaSpec = AnalysisSpecification(true, true, true, true, true, true, false, false)
testAnaSpecSym = AnalysisSpecification(true, true, true, true, true, true, false, true)



@testset "BellDiagonalQudits.jl" begin

    @test length(createRandomCoordStates(10, 3, "magicSimplex")) == 10
    @test length(createStandardIndexBasis(4, 10).basis) == 16
    @test testDictionaries3[1][2, 1] == 6
    @test createDensityState(CoordState(1 / 9 * ones(9), "UNKNOWN"), testStandardBasis3).coords ≈ 1 / 9 * ones(9)
    @test length(createBipartiteWeylOpereatorBasis(3)) == 81
    @test length(vertices_list(tovrep(testStandardKernel3))) == 12
    @test length(generateAllSymmetries(testStandardBasis3, 3)) == 432
    @test getSymCoords(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        generateAllSymmetries(testStandardBasis3, 3)[1:2]
    ) ==
          [
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    ]
    @test analyseCoordState(
        3,
        CoordState(1 / 9 * ones(9), "UNKNOWN"),
        testAnaSpec,
        testStandardBasis3,
        testStandardKernel3,
        testBipartiteWeylBasis3,
        testDictionaries3,
        testMub3,
        missing
    ) == AnalysedCoordState(
        CoordState(1 / 9 * ones(9), "UNKNOWN"),
        true,
        true,
        true,
        false,
        false,
        false,
        missing
    )

end
