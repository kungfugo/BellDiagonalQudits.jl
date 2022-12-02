[![CI][ci-img]][ci-url]
[![codecov][codecov-img]](codecov-url)
[![][docs-dev-img]][docs-dev-url]

[ci-img]: https://kungfugo.github.io/BellDiagonalQudits.jl/actions/workflows/CI.yml/badge.svg
[ci-url]: https://kungfugo.github.io/BellDiagonalQudits.jl/actions/workflows/CI.yml
[cov-img]: http://codecov.io/github/kungfugo/BellDiagonalQudits.jl/coverage.svg
[cov-url]: https://codecov.io/github/kungfugo/BellDiagonalQudits.jl
[docs-dev-url]: https://kungfugo.github.io/BellDiagonalQudits.jl/dev/

# BellDiagonalQudits.jl

_Generate and analyze Bell diagonal Qudits with Julia_

A package for generation and entanglement classification of Bell diagonal quantum states.

Bell diagonal states are generated as mixtures maximally entangled Bell states, which are related by Weyl transformations. The special propterties of these states, e.g. symmetries, allow efficient methods to be leveraged for the detection of entanglement, including its generally hard to detect form of PPT/bound entanglement.

This package provides methods to sample states, to numerically generate entanglement witnesses and to apply and extend further criteria to detect entanglement or separability in general dimension. For a precise description of implemented methods and related research results see [1], [2], [3] and the references therein.

## Package Features

- Create mixtures of maximally entangled Bell states based on Weyl transformations in any dimension
- Classify Bell diagonal states as separable, PPT/bound entangled or NPT/free entangled
- Generate numerical entanglement witnesses for Bell diagonal states
- Generate entanglement conserving symmetries and use them for entanglement classification

## Installation

BellDiagonalQudits can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```
pkg> add BellDiagonalQudits
```

The package can be loaded via

```julia
julia> using BellDiagonalQudits
```

## Documentation

Documentation is available at [https://kungfugo.github.io/BellDiagonalQudits.jl/dev/](https://kungfugo.github.io/BellDiagonalQudits.jl/dev/)

## References

[1] Popp, C., Hiesmayr, B.C., _Almost complete solution for the NP-hard separability problem of Bell diagonal qutrits_, Sci Rep 12, 12472 (2022), [https://doi.org/10.1038/s41598-022-16225-z](https://doi.org/10.1038/s41598-022-16225-z)

[2] Baumgartner, B., Hiesmayr, B.C., Narrenhofer, H. _A special simplex in the state space for entangled qudits_, J. Phys. A Math. Theor. 40, 7919 (2007), [https://doi.org/10.1088/1751-8113/40/28/s03] (https://doi.org/10.1088/1751-8113/40/28/s03)

[3] Popp, C., Hiesmayr, B.C., _Bound Entanglement of Bell Diagonal Pairs of Qutrits and Ququarts: A Comparison_, arXiv (2022), [https://arxiv.org/abs/2209.15267] (https://arxiv.org/abs/2209.15267)
