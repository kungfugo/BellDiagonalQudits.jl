# Manual

This manual shows how to use the package by using it to sample a set of uniformly distributed Bell diagonal states and to analyse their entanglement properties. Here, we apply criteria for separability and entanglement to determine the entanglement class of generated bipartite qutrits, i.e. `d=3`. In this system, bound entangled states, i.e. entangled states with positive partial transposition (PPT) that cannot be used for entanglement distillation, exist. The entanglement classes are labeled "SEP" for separability, "BOUND" for bound entanglement, "FREE" for entanglement with negative partial transposition (i.e. distillable) and "PPT_UNKNOWN" for PPT states that could not be classified as entangled or separable.

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

Create a basis of maximally entangled bipartite Bell states in `d^2` dimensions. Each Bell basis state is created by applying a certain Weyl transformation to the maximally entangled state. Sample Bell diagonal random mixed states of those Bell states, represented by their `d^2` coordinates (mixing probabilities) in the Bell basis.

\
**Bell basis generation**

Create a Bell basis `myBasis` by applying each of the `d^2` Weyl transformations `W_{k,l} \\otimes \\mathbb{1}_d` to the maximally entangled state.`myBasis.basis` contains the enumerated Bell states in computational basis together with the indices of the corresponding Weyl transformation. `myBasisDict` contains the dictionaries to relate the enumerated `d^2` Bell basis states to the double indices `(k,l)` of the corresponding Weyl transformation.

```julia
d = 3
myBasis = create_standard_indexbasis(d,10)
myBasisDict = create_dictionary_from_basis(myBasis)
```

\
**State sampling**

Create uniformly distributed random representations of quantum states by specifying the coordinates of the state in the created Bell basis. The coordinates represent the mixing probabilities of the Bell basis states. Here, we create only states with Bell coordinates within the "enclosure polytope", which is defined by the limitation of all coordinates (mixing probabilities) to be smaller than or equal to `1/d`. This subset is known to contain all (but not only) states with positive partial transposition (PPT), which can be separable or bound entangled.

```julia
myCoordStates = uniform_bell_sampler(100, d, :enclosurePolytope)
```

Create `DensityState`s including the density matrix of each state in the computational basis, created by mixing the Bell states of `myBasis` according to the `coords` of `myCoordStates`.

```julia
myDensityStates = map(x->create_densitystate(x, myBasis), myCoordStates)
```

## Analysis prerequisites

Now create the analysis objects required for the entanglement classification using several criteria for entanglement or separability.

\
**Separable kernel polytope**

The kernel polytope is known to contain only Bell coordinates that represent separable states. It is defined as the convex hull of vertices related to special separable states called "subgroup states". The related `kernel check` tests if the Bell coordinates of a given unclassified state are contained in this convex hull and thus indicates separability.

```julia
mySepKernel = create_kernel_polytope(d, myBasis)

```

If additional separable states `newSepDensityStates` are known, the kernel polytope can be exteded to a larger convex hull in order to improve the kernel check for separability.
For a (trivial) example, consider the separable, maximally mixed state having all Bell states mixed equally with probability `1/d^2`. First,specify the coordinates in the Bell basis and set the `eClass` of the corresponding `CoordState` to "SEP". Then, calculate the density matrix and create the `DensityState`. Finally, extend the kernel polytope `mySepKernel` by the array containing this separable state.

```julia
maxMixedCoordState = CoordState(1/d^2*ones(d^2), "SEP")
maxMixedDensityState = create_densitystate(maxMixedCoordState, myBasis)
newSepDensityStates = [maxMixedDensityState]
myExtendedKernel = extend_vpolytope_by_densitystates(tovrep(mySepKernel), newSepDensityStates, 10)

```

\
**Weyl operator basis**

Use the Weyl operators to construct a basis of the space of `(d^2,d^2)` matrices. This object is used for the `spinrep check` indicating separability according to the representation of the density matrix of a given state in this basis.

```julia
myWeylOperatorBasis = create_bipartite_weyloperator_basis(d)
```

\
**Mutually unbiased bases (MUBs)**

Create the a set of mutually unbiased bases (MUBs) constructed with the Weyl operators and represented in the computational basis.

```julia
myMub = create_standard_mub(d)
```

\
**Symmetries**

Generate entanglement class conserving symmetries represented as permutations of state coordinates in the Bell basis. Given a classified state, the orbit, i.e. the set of states that are generated by applying all symmetries to the classified state, is known to be of the same entanglement class. The symmetries can be used to improve the entanglement classification.

```julia
mySyms = generate_symmetries(myBasis, d)
```

\
**Entanglement witnesses**

Generate `n` numerical entanglement witnesses by numerical optimization over the set of separable states. Here, the entanglement witnesses are represented by their coordinates in the Bell basis, an upper, and a lower bound. For all separable states, the inner product of the state and witness coordinates obeys these bounds. A violation of the inner product of an unknown state and the witness thus indicates entanglement. Use `iterations` runs to improve the determined upper and lower bounds. Other optimization methods than the default `NelderMead` can be used.

```julia
n = 2
myOptimizedEWs = create_random_bounded_ews(
    d,
    myBasis,
    n,
    true,
    20
    )
```

```julia
myOptimizedCoodEWs = map(x->get_bounded_coordew(x), myOptimizedEWs)
```

## Entanglement classification

**Analysis specification**

Specify, which entanglement checks to use. See properties of type `AnalysisSpecification`. In this case we check separability with the kernel and spinrep check and test for entanglement using the ppt, realignment, concurrence_qp and numeric_ew check.

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

If `useSymmetries == false` in the analysis specification `myAnaSpec` use `analyse_coordstate`, else use `sym_analyse_coordstate` to leverage the symmetries `mySyms` for improved classification.

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

## Entanglement distillation

Create and distill a Bell-diagonal state with the FIMAX protocol. First, create a test state.

```julia
d=3
testBDS = create_densitystate(CoordState([0.5, 0.5 / 8, 0.5 / 8, 0.5 / 8, 0.5 / 8, 0.5 / 8, 0.5 / 8, 0.5 / 8, 0.5 / 8], "UNKNOWN"), myBasis).densityMatrix

```

To execute one iteration of the FIMAX routine, run:

```julia
FIMAX_routine_results = FIMAX_routine(testBDS, 2, d, myBasis)
```

To iterate this procedure until a target fidelity of 0.99 with the maximally entangled state is achieved, execute:

```julia
FIMAX_protocol_results = iterative_FIMAX_protocol(testBDS, 0.99, 2, d, myBasis, 100)
```