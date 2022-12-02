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

Create a basis of maximally entangled Bell states in `d` dimensions, indexed by the corresponding Weyl operator and sample some random states represented by `d^2` coordinates in the "enclosurePolytope", which is known to contain all PPT states.

\
**Bell basis generation**

Create a basis and a dictionary to relate the `d^2` coordinates of a Bell diagonal state to the double indices `(k,l)` of the corresponding Weyl operator `W_{k,l}`.

```julia
myBasis = create_standard_indexbasis(d,10)
myBasisDict = create_dictionary_from_basis(myBasis)
```

\
**State sampling**

Create uniformly distributed random representations of quantum states by specifying the coordinates of the state in the created Bell basis.

```julia
myCoordStates = uniform_bell_sampler(100, d, "enclosurePolytope")
```

Create `DensityState`s including the density matrix represented in the computational basis.

```julia
myDensityStates = map(x->create_densitystate(x, myBasis), myCoordStates)
```

## Analysis prerequisites

Now create the analysis objects required for the entanglement classification.

\
**Separable kernel polytope**

The kernel polytope is known to contain only Bell coordinates that relate to separable states.

```julia
mySepKernel = create_kernel_polytope(d, myBasis)
```

If additional separable states `newSepDensityStates` are known, the kernel polytope can be exteded to improve the kernel check for separability.

```julia
myExtendedKernel = extend_vpolytope_by_densitystates(tovrep(mySepKernel), newSepDensityStates)
```

\
**Weyl operator basis**

Use the Weyl operators to construct a basis of the space of `(d^2,d^2)` matrices.

```julia
myWeylOperatorBasis = create_bipartite_weyloperator_basis(d)
```

\
**Mutually unbiased bases (MUBs)**

Create the MUBS using the Weyl operators to construct them from the computational basis.

```julia
myMub = create_standard_mub(d)
```

\
**Symmetries**

Generate entanglement class conserving symmetries represented as permutations of state coordinates in the Bell basis.

```julia
mySyms = generate_symmetries(myBasis, d)
```

\
**Entanglement witnesses**

Generate `n` numerical entanglement witnesses by numerical optimization over the set of separable states. Use `iterations` runs to improve the determined upper and lower bounds. Other optimization methods than the default `NelderMead` can be used.

```julia
myOptimizedEWs = create_random_bounded_ews(
    d,
    myBasis,
    n,
    true,
    50
    )
```

Represent entanglement witnesses via their coordinates in Bell Basis.

```julia
myOptimizedCoodEWs = map(x->get_bounded_coordew(x), myOptimizedEWs)
```

## Entanglement classification

**Analysis specification**

Specify, which entanglement checks to use. See properties of type `AnalysisSpecification`.

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

If `useSymmetries == false` in the analysis specification `myAnaSpec` use `analyse_coordstate`, else use `sym_analyse_coordstate` to include the symmetry analysis.

```julia
f(x) = analyse_coordstate(
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
classify_analyzed_states!(myAnalysedCoordStates)
```

Identify e.g. bound entangled states as

```julia
myBoundStates = filter(x->x.coordState.eClass == "BOUND", myAnalysedCoordStates)
```
