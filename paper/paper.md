--- 
title: 'BellDiagonalQudits: A package for entanglement analyses of mixed
maximally entangled qudits' tags:

- Julia
- quantum information
- quantum physics
- entanglement
- separability problem
- Bell states
  authors:
- name: Christopher Popp orcid: 0000-0002-2352-3731 equal-contrib: true
  affiliation: 1
  affiliations:
- name: Faculty of Physics, University of Vienna, Währingerstrasse 17, 1090,
  Vienna, Austria index: 1

date: 20 September 2022 
bibliography: paper.bib 
---

# Summary

In the field of quantum information and technology, entanglement of quantum
states called qudits is regarded as ressource for quantum and classical
information processing tasks and allows the use of algorithms with better
performance than any classical algorithm for certain applications like
superdense coding, teleportation or computing. Despite its theoretical and
practial relevance, there is now general method to determine whether a given
quantum state is entangled or not due to a special and hard to detect form of
entanglement called bound entanglement. Bipartite Bell states are maximally
entangled states of two qudits and are of high relevance for application in
quantum technologies. Mixtures of those states can be entangled, including it's
bound form, or separable. Leveraging their special properties, Bell states can
be numerically generated and analyzed for their entanglement properties by
various methods implemented in this Julia package.

# Statement of need

`BellDiagonalQudits` is an Quantum Information-affiliated Julia package for generation and
entanglement analyses of mixed maximally entangled bipartite qudits in general
dimension. The API for `BellDiagonalQudits` provides an user-friendly interface
to generate representations of Bell diagonal quantum states and to analyze their
entanglement properties with various general or specialized criteria to detect
entanglement or separability. Leveraging geometric properties of a certain class
of mixed Bell states that are related by Weyl transformations,
`BellDiagonalQudits` combines known analytical results [@baumgartner] and
numerical methods for quantum state representation and analysis. It leverages
and depends on the Julia packages `QuantumInformation`, convex sets of
`LazySets` [@lazysets] and optimization of `Optim` [@optim].

`BellDiagonalQudits` was designed to be used by researchers in quantum science
and quantum information theory. It has already been used in multiple scientific
publications [@PoppACS; @PoppBoundEntComparison] in the context of entanglement
classification and detection of bipartite qudits in dimension three and
four. The combination of efficient state generation via random sampling or
determenistic procedures and implementation of both frequently used and
specialized entanglement and separability detectors supports the research of
entanglement in several ways. From a general point of view, entangled Bell
states are well accessible for powerfull methods of entanglement and
separability detection, leveraging their symmetries and geometric properties. It
was shown in [@PoppACS; @PoppBoundEntComparison; @hiesmayr] that a significant
share of the analyzed group of mixed Bell states related by Weyl transformations
are bound entangled, offering a systematic way to generate and investigatie
those states with respect to the separability problem in different
dimensions. In addition to general methods applicable to any Bell diagonal states,
`BellDiagonalQudits` provides features to generate the special symmetries of Bell
states that are related via Weyl transformations. These symmetries are leveraged
for entanglement classification and the numerical generation of specialized
entanglement witnesses in any dimension. Furthermore, the implemented methods
of `BellDiagonalQudits` can be used in various quantum information processing
tasks involving Bell diagonal states in any dimension like Quantum Key Distribution
or entanglement verification. `BellDiagonalQudits` uses and integrates well with
the general interface of `QuantumInformation` allowing the investigation of Bell
states in the context of quantum channels, entanglement measures or entropy.

# Relation to research projects

The methods of `BellDiagonalQudits` to generate and analyze Bell diagonal states
in general dimension are based on analytical properties summarized in
[@baumgartner]. Extensions of those methods and efficient implementation in
`BellDiagonaQudits` enabled the detailed analysis [@PoppACS] of Bell diagonal
qudits in three dimensions (qutrits), providing an operational solution to the
separability problem for those states. Additionally the relative shares of
separable and (bound) entangled states were precisely determined among the Bell
diagonal states. In [@PoppBoundEntComparison], higher dimensions were
considered, focusing on a detailed comparison and geometric properties of
separable states in dimension three and four.

# Package information

`BellDiagonalQudits` is available on Github at [Link
repository](https://github.com/kungfugo/BellDiagonalQudits.jl). The package
documentation is available at [Link
documentation](https://kungfugo.github.io/BellDiagonalQudits.jl/dev/) and
provides examples of usage.

# Acknowledgements

I acknowledge support from Beatrix C. Hiesmayr for review and validation of
implemented methods.

# References
