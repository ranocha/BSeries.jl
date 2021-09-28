# BSeries.jl

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
@misc{ranocha2021bseries,
  title={{BSeries.jl}: {C}omputing with {B}-series in {J}ulia},
  author={Ranocha, Hendrik and Ketcheson, David I},
  year={2021},
  month={09},
  howpublished={\url{https://github.com/ranocha/BSeries.jl}},
  doi={}
}
```
Please also cite the appropriate references for specific functions you use,
which can be obtained from their docstrings.


## License and contributing

This project is licensed under the MIT license (see [License](@ref)).
Since it is an open-source project, we are very happy to accept contributions
from the community. Please refer to the section [Contributing](@ref) for more
details.
