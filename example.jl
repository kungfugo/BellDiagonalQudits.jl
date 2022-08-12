include("src/MagicSimplex.jl")

using .MagicSimplex
using JLD
using LazySets

# Load exemplary data for analysis 
myKernel = load("Data/standardKernelVPt3.jld", "sepHull");
myEWs = load("Data/boundedCoordEWs.jld", "boundedCoordEWs");
myWeylOpBasis = load("Data/bipartiteWeylOpBasis3.jld", "bipartiteWeylOpBasis");

# Create basis 
myBasis = createStandardIndexBasis(3, 10);
myBasisDict = MagicSimplex.createDictionaryFromBasis(myBasis);

# Symmetries 
mySyms = generateAllSymmetries(myBasis, 3);

# Load MUBs 
myMubs = MagicSimplex.createStandardMub(3);

# Create random states
myCoordStates = createRandomCoordStates(100, 3, "enclosurePolytope");

# Specify what to analyse
myAnaSpec = AnalysisSpecification(
    true, true, true, true, true, true, true, true
)

#Run analysis 
if myAnaSpec.useSymmetries
    f(x) = symAnalyseCoordState(
        3,
        x,
        mySyms,
        myAnaSpec,
        myBasis,
        myKernel,
        myWeylOpBasis,
        myBasisDict,
        myMubs,
        myEWs
    )
else
    f(x) = analyseCoordState(
        3,
        x,
        myAnaSpec,
        myBasis,
        myKernel,
        myWeylOpBasis,
        myBasisDict,
        myMubs,
        myEWs
    )
end

myAnalysedStates = map(x -> f(x), myCoordStates)

# Classify analysed states 
classifyAnalyzedStates!(myAnalysedStates)