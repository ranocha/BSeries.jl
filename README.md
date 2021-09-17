# BSeries

[![Build Status](https://github.com/ranocha/BSeries.jl/workflows/CI/badge.svg)](https://github.com/ranocha/BSeries.jl/actions?query=workflow%3ACI)
[![Coverage Status](https://coveralls.io/repos/github/ranocha/BSeries.jl/badge.svg?branch=main)](https://coveralls.io/github/ranocha/BSeries.jl?branch=main)
[![codecov](https://codecov.io/gh/ranocha/BSeries.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ranocha/BSeries.jl)

A collection of functionality around [B-series](https://en.wikipedia.org/wiki/Butcher_group) in [Julia](https://julialang.org/). See

- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)


## API Documentation

Please note that this project is work in progress.

BSeries.jl re-exports everything from
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl).
However, if you rely on functionality from that package,
you should also include it explicitly in your project dependencies
to track breaking changes, since the version numbers of RootedTrees.jl
and BSeries.jl are not necessarily synchronized.
