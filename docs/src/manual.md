# Manual

## Package installation
BellDiagonalQudits can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run 

```  
pkg> add BellDiagonalQudits
```
The package can be loaded via
```julia   
julia> using BellDiagonalQudits
```

## State generation
First, create a basis of maximally entangled Bell states in `d` dimensions, indexed by the corresponding Weyl operator and sample some random states represented by ``d^2`` coordinates in the "enclosurePolytope", which is known to contain all PPT states.

\
**Bell basis generation**

Create a basis and a dictionary to relate the ``d^2`` coordinates of a Bell diagonal state to the double indices ``(k,l)`` of the corresponding Weyl operator ``W_{k,l}``. 

```julia   
myBasis = createStandardIndexBasis(d)
myBasisDict = createDictionaryFromBasis(myBasis)
```

\
**State sampling**

Create representations of quantum states by specifying the coordinates of the state in the created Bell basis.

```julia   
myCoordStates = createRandomCoordStates(100, d, "enclosurePolytope")
```

Create `DensityState`s including the density matrix represented in the computational basis.

```julia   
myDensityStates = map(x->createDensityState(x, myBasis), myCoordStates)
```

## Analysis prerequisites
Now create the analysis objects required for the entanglement classification.

\
**Separable kernel polytope**

The kernel polytope is known to contain only points that relate separable states via their Bell coordinates.

```julia       
mySepKernel = createKernelPolytope(d, myBasis)
```

If additional separable `DensityState`s `newSepDensityStates` are known, the kernel polytope can be exteded to improve the kernel check for separability.

```julia   
myExtendedKernel = extendSepVPolytopeBySepStates(tovrep(mySepKernel), newSepDensityStates)
```

\
**Weyl operator basis**

Use the Weyl operators to construct a basis of the space of `(d^2,d^2)` matrices.

```julia   
myWeylOperatorBasis = createBipartiteWeylOpereatorBasis(d)    
```
\
**Mutually unbiased bases (MUBs)**

Create the standard MUBS generated via the Weyl operators.

```julia   
myMub = createStandardMub(d)
```

\
**Symmetries**

Generate entanglement class conserving symmetries represented as permutations of state coordinates in the Bell basis.

```julia   
mySyms = generateAllSymmetries(myBasis)
```

\
**Entanglement witnesses**

Generate `n` numerical entanglement witnesses by numerical optimization over the set of separable states. Use `iterations` runs to improve the determined upper and lower bounds. Other optimization methods than `NelderMead` can be used.

```julia   
myOptimizedEWs = createRandomBoundedWits(
    d,
    myBasis,
    n,
    true,
    50,
    NelderMead
)
```
Represent entanglement witnesses via their coordinates in Bell Basis.

```julia   
myOptimizedCoodEWs = map(x->getBoundedCoordEw(x), myOptimizedEWs)
```
## Entanglement classification 

**Analysis specification**

Specify which entanglement checks to use. See properties of type `AnalysisSpecification`.
 ```julia   
myAnaSpec = AnalysisSpecification(
    true,
    true,
    true,
    true,
    true,
    false,
    true,
    false
)
```

\
**Apply analysis to all generated states**

If `useSymmetries == false` in the analysis specification `myAnaSpec` use `analyseCoordState`, else use `symAnalyseCoordState`.
```julia   
f(x) = analyseCoordState(
    d,
    x,
    myAnaSpec,
    myBasis,
    mySepKernel,
    myWeylOperatorBasis,
    myBasisDict,
    missing,
    myOptimizedCoodEWs
)

    myAnalysedCoordStates = map(x->f(x), myCoordStates)
```

Finally use analysis results to set `CoordState.eClass` to assign the entanglement class to the states.

```julia
classifyAnalyzedStates!(myAnalysedStates)
```