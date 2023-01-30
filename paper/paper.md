---
title: "BellDiagonalQudits: A package for entanglement analyses of mixed maximally entangled qudits"
tags:
  - Julia
  - quantum information
  - quantum physics
  - entanglement
  - separability problem
  - Bell states
authors:
  - name: Christopher Popp
    orcid: 0000-0002-2352-3731
    affiliation: 1
affiliations:
  - name: Faculty of Physics, University of Vienna, Währingerstrasse 17, 1090, Vienna Austria
    index: 1
date: 15 December 2022
bibliography: paper.bib
---

# Summary

In the field of quantum information and technology, entanglement in quantum
systems called qudits is regarded as resource for quantum-based or quantum-assisted
information processing tasks. It allows new ways of information processing like
quantum computation or quantum teleportation and provides possibilities for
speedups in a variety of algorithms with applications like search or optimization.
Despite its theoretical and practical relevance, there is no general method to determine whether a given
quantum state is entangled or not due to a special and hard to detect form of
entanglement called bound entanglement. Bipartite Bell states are maximally
entangled states of two qudits and are of high relevance for application in
quantum technologies. Mixtures of those states can be entangled, including its
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
`BellDiagonalQudits` combines known analytical results by @baumgartner and
numerical methods for quantum state representation and analysis. It leverages
and depends on the Julia package `QuantumInformation`, the convex sets of
`LazySets` [@lazysets] and the optimization methods of `Optim` [@optim].

`BellDiagonalQudits` was designed to be used by researchers in quantum science
and quantum information theory. It has already been used in multiple scientific
publications, e.g. in @PoppACS [@PoppACS:2022] and @PoppBoundEntComparison:2022 in the context of entanglement
classification and detection of bipartite qudits in dimension three and
four. The combination of efficient state generation via random sampling or
deterministic procedures and implementation of both frequently used and
specialized entanglement and separability detectors supports the research of
entanglement in several ways. From a general point of view, entangled Bell
states are well accessible for powerful methods of entanglement and
separability detection, leveraging their symmetries and geometric properties. It
was shown in that a significant share of the group of mixed Bell states related by Weyl transformations
are bound entangled [@PoppACS;@PoppBoundEntComparison;@hiesmayr], offering a systematic way to generate and investigate
those states with respect to the separability problem in different
dimensions. In addition to general methods applicable to any Bell diagonal state,
`BellDiagonalQudits` provides features to generate the special symmetries of Bell
states that are related via Weyl transformations. These symmetries are leveraged
for improved entanglement classification and the numerical generation of specialized
entanglement witnesses in any dimension. Furthermore, the implemented methods
of `BellDiagonalQudits` can be used in various quantum information processing
tasks involving Bell diagonal states in any dimension like Quantum Key Distribution
or entanglement verification. `BellDiagonalQudits` uses and integrates well with
the general interface of `QuantumInformation` allowing the investigation of Bell
states in the context of quantum channels, entanglement measures or entropy.

# Relation to research projects

The methods of `BellDiagonalQudits` to generate and analyze Bell diagonal states
in general dimension are based on analytical properties summarized by
@baumgartner:2007. Extensions of those methods and efficient implementation in
`BellDiagonalQudits` enabled the detailed analysis by @PoppACS:2022 of Bell diagonal
qudits in three dimensions (qutrits), providing an operational solution to the
separability problem for those states. Additionally the relative shares of
separable and (bound) entangled states were precisely determined among the Bell
diagonal states. In @PoppBoundEntComparison:2022, higher dimensions were
considered, focusing on a detailed comparison and geometric properties of
separable states in dimension three and four.

# Package information

`BellDiagonalQudits` is available on Github at [https://github.com/kungfugo/BellDiagonalQudits.jl](https://github.com/kungfugo/BellDiagonalQudits.jl). The package
documentation is available at [https://kungfugo.github.io/BellDiagonalQudits.jl/dev/](https://kungfugo.github.io/BellDiagonalQudits.jl/dev/) and
provides examples of usage.

# Acknowledgments

I acknowledge support from Beatrix C. Hiesmayr for review and validation of
implemented methods.

# References
