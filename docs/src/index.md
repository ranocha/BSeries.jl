# BSeries.jl

[BSeries.jl](https://github.com/ranocha/BSeries.jl)
is a collection of functionality around
[B-series](https://en.wikipedia.org/wiki/Butcher_group)
in [Julia](https://julialang.org/). See for example

- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)

BSeries.jl re-exports everything from
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl).
However, if you rely on functionality from that package,
you should also include it explicitly in your project dependencies
to track breaking changes, since the version numbers of RootedTrees.jl
and BSeries.jl are not necessarily synchronized.

The main API of BSeries.jl consists of the following components.

- B-series behave like `AbstractDict`s mapping
  (abstract) `RootedTree`s to coefficients.
- The B-series of time integration methods such as Runge-Kutta methods
  can be constructed by the function [`bseries`](@ref).
- The algebraic structures of the composition law and the substitution law are
  implemented via [`compose`](@ref) and [`substitute`](@ref).
- Backward error analysis can be performed via
  [`modified_equation`](@ref)s and [`modifying_integrator`](@ref)s.

Further information is provided in the following tutorials and
API documentation.


## Installation

[BSeries.jl](https://github.com/ranocha/BSeries.jl)
is a registered Julia package. Thus, you can install it from the Julia REPL via
```julia
julia> using Pkg; Pkg.add("BSeries")
```

If you want to update BSeries.jl, you can use
```julia
julia> using Pkg; Pkg.update("BSeries")
```
As usual, if you want to update BSeries.jl and all other
packages in your current project, you can execute
```julia
julia> using Pkg; Pkg.update()
```


## Referencing

If you use
[BSeries.jl](https://github.com/ranocha/BSeries.jl)
for your research, please cite it using the bibtex entry
```bibtex
@online{ketcheson2021computing,
  title={Computing with {B}-series},
  author={Ketcheson, David I and Ranocha, Hendrik},
  year={2021},
  month={11},
  eprint={2111.11680},
  eprinttype={arXiv},
  eprintclass={math.NA}
}
```
In addition, you can also refer to BSeries.jl directly as
```bibtex
@misc{ranocha2021bseries,
  title={{BSeries.jl}: {C}omputing with {B}-series in {J}ulia},
  author={Ranocha, Hendrik and Ketcheson, David I},
  year={2021},
  month={09},
  howpublished={\url{https://github.com/ranocha/BSeries.jl}},
  doi={10.5281/zenodo.5534602}
}
```
Please also cite the appropriate references for specific functions you use,
which can be obtained from their docstrings.


## License and contributing

This project is licensed under the MIT license (see [License](@ref)).
Since it is an open-source project, we are very happy to accept contributions
from the community. Please refer to the section [Contributing](@ref) for more
details.
