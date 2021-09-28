# BSeries

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ranocha.github.io/BSeries.jl/stable)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ranocha.github.io/BSeries.jl/dev)
[![Build Status](https://github.com/ranocha/BSeries.jl/workflows/CI/badge.svg)](https://github.com/ranocha/BSeries.jl/actions?query=workflow%3ACI)
[![Coveralls](https://coveralls.io/repos/github/ranocha/BSeries.jl/badge.svg?branch=main)](https://coveralls.io/github/ranocha/BSeries.jl?branch=main)
[![Codecov](https://codecov.io/gh/ranocha/BSeries.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ranocha/BSeries.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
<!-- [![DOI](https://zenodo.org/badge/DOI/.svg)](https://doi.org/) -->

A collection of functionality around [B-series](https://en.wikipedia.org/wiki/Butcher_group) in [Julia](https://julialang.org/). See

- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)


## API Documentation

The Julia library
[BSeries.jl](https://github.com/ranocha/BSeries.jl)
is work in progress. Nevertheless, we follow semantic versioning. Thus, you can
safely use the package right now. Extended documentation will be provided in the
future.

BSeries.jl re-exports everything from
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl).
However, if you rely on functionality from that package,
you should also include it explicitly in your project dependencies
to track breaking changes, since the version numbers of RootedTrees.jl
and BSeries.jl are not necessarily synchronized.
