# Code generation
This tutorial shows how to generate C code to compute expressions found using [BSeries.jl](https://github.com/ranocha/BSeries.jl).  Although [BSeries.jl](https://github.com/ranocha/BSeries.jl) is compatible with three [symbolic backends](https://ranocha.de/BSeries.jl/dev/tutorials/symbolic_computations/), it's currently easiest to perform code generation using [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

First, we generate the B-series that we want work with.  Here we take a generic 2nd-order RK method and generate terms only up to 3rd order, in order to work with the leading truncation error.


```julia
using BSeries, Latexify, Symbolics

@variables α
A = [0 0; 1/(2*α) 0]; b = [1-α, α]; c = [0, 1/(2*α)]

RK22 = bseries(A, b, c, 3)
exact = ExactSolution(RK22)
truncation_error = RK22-exact
```

Next we set up the ODE of interest, and evaluate the B-series with that right-hand side.


```julia
@variables h
@variables p q
f = [p * (2 - q), q * (p - 1)]

du = evaluate(f,[p,q],h,RK22-exact)
```

Finally, we generate a C function that evaluates the expressions above.


```julia
build_function(du, α, p, q, h, target=Symbolics.CTarget())
```

The [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) routine `build_function` can also generate code in Julia, MATLAB, and Stan; see the [documentation](https://symbolics.juliasymbolics.org/stable/manual/build_function/#build_function) for details and other options.
