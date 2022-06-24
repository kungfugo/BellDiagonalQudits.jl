using MagicSimplex

# Load exemplary data for analysis 
myKernel = load("MagicSimplex/Data/standardKernelVPt3.jld", "sepHull");
myEWs = load("MagicSimplex/Data/boundedCoordEWs.jld", "boundedCoordEWs");
myWeylOpBasis = load("MagicSimplex/Data/bipartiteWeylOpBasis3.jld", "bipartiteWeylOpBasis");

# Create basis 
myBasis = createStandardIndexBasis(3, 10);
myBasisDict = MagicSimplex.createDictionaryFromBaxsis(myBasis);

# Create random states
myCoordStates = createRandomCoordStates(100, 3, "enclosurePolytope");

# Symmetries 
mySyms = generateAllSymmetries(myBasis, 3)

# Load MUBs 
myMubs = MagicSimplex.createStandardMub(3);

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

