# [Benchmarks](@id benchmarks)

Here, we collect some simple benchmarks of
[BSeries.jl](https://github.com/ranocha/BSeries.jl).
Take them with a grain of salt since they run on virtual machines in the
cloud to generate the documentation automatically. You can of course also
copy the code and run the benchmarks locally yourself.


## [Comparing different symbolic packages](@id benchmarks-symbolic-computations)

Symbolic computations of [`modified_equation`](@ref)s and
[`modifying_integrator`](@ref)s in
[BSeries.jl](https://github.com/ranocha/BSeries.jl)
support

- [SymEngine.jl](https://github.com/symengine/SymEngine.jl),
- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), and
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

as symbolic backends. Here, we compare them in the context of the explicit
midpoint method and the nonlinear oscillator ODE

```math
u'(t) = \frac{1}{\| u(t) \|^2} \begin{pmatrix} -u_2(t) \\ u_1(t) \end{pmatrix}.
```

This particular combination of explicit Runge-Kutta method and ODE is special
since the explicit midpoint method is unconditionally energy-conserving for this
problem[^RanochaKetcheson2020].

First, we set up some code to perform the benchmarks. Here, we use a very naive
approach, run the code twice (to see the effect of compilation) and use `@time`
to print the runtime. More sophisticated approaches should of course use
something like `@benchmark` from
[BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl).
However, this simple and cheap version suffices to compare the orders of
magnitude.

```@example benchmark-nonlinear-oscillator
using BSeries, StaticArrays

function benchmark(u, dt, subs, order)
  # explicit midpoint method
  A = @SArray [0 0; 1//2 0]
  b = @SArray [0, 1//1]
  c = @SArray [0, 1//2]

  # nonlinear oscillator
  f = [-u[2], u[1]] / (u[1]^2 + u[2]^2)

  println("\n Computing the series coefficients:")
  @time coefficients = modifying_integrator(A, b, c, order)
  @time coefficients = modifying_integrator(A, b, c, order)

  println("\n Computing the series including elementary differentials:")
  @time series = modifying_integrator(f, u, dt, A, b, c, order)
  @time series = modifying_integrator(f, u, dt, A, b, c, order)

  substitution_variables = Dict(u[1] => 1//1, u[2] => 0//1)

  println("\n Substituting the initial condition:")
  @time subs.(series, (substitution_variables, ))
  @time subs.(series, (substitution_variables, ))

  println("\n")
end
```

Next, we load the symbolic packages and run the benchmarks.

```@setup benchmark-nonlinear-oscillator
using SymPy # generates annoying output online when conda installs sympy
```

```@example benchmark-nonlinear-oscillator
using SymEngine: SymEngine
using SymPy: SymPy
using Symbolics: Symbolics

println("SymEngine")
dt   = SymEngine.symbols("dt")
u    = SymEngine.symbols("u1, u2")
subs = SymEngine.subs
benchmark(u, dt, subs, 8)

println("SymPy")
dt   = SymPy.symbols("dt")
u    = SymPy.symbols("u1, u2")
subs = SymPy.subs
benchmark(u, dt, subs, 8)

println("Symbolics")
Symbolics.@variables dt
u = Symbolics.@variables u1 u2
subs = Symbolics.substitute
benchmark(u, dt, subs, 8)
```

These results were obtained using the following versions.

```@example
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(["BSeries", "RootedTrees", "SymEngine", "SymPy", "Symbolics"],
           mode=PKGMODE_MANIFEST)
nothing # hide
```


## [Comparison with other packages](@id benchmarks-other-packages)

There are also other open source packages for B-series. Currently, we are aware
of the Python packages

- [`BSeries`](https://github.com/ketch/BSeries)
- [`pybs`](https://github.com/henriksu/pybs)

If you know about similar open source packages out there, please inform us, e.g.,
by [creating an issue](https://github.com/ranocha/BSeries.jl/issues/new/choose)
on GitHub.

The packages listed above and [BSeries.jl](https://github.com/ranocha/BSeries.jl)
all use different approaches and have different features. Thus, comparisons
must be restricted to their common subset of features. Here, we present some
simple performance comparisons. Again, we just use (the equivalent of) `@time`
twice to get an idea of the performance after compilation, allowing us to
compare orders of magnitude.


First, we start with the Python package
[`BSeries`](https://github.com/ketch/BSeries)
and the following benchmark script.

````@example
using BSeries, Markdown                                                   # hide
filename = joinpath(pathof(BSeries) |> dirname |> dirname, "docs", "src", # hide
  "benchmark_python_bseries.py")                                          # hide
script = "```python\n"                                                    # hide
for line in Iterators.drop(readlines(filename), 4)                        # hide
  startswith(line, "with") && continue                                    # hide
  line = replace(line, "  print" => "print")                              # hide
  global script                                                           # hide
  script = script * replace(line, ", file=io" => "") * "\n"               # hide
end                                                                       # hide
script = script * "```\n"                                                 # hide
Markdown.parse(script)                                                    # hide
````

The results are as follows.

````@example
using BSeries, Markdown                                                   # hide
filename = joinpath(pathof(BSeries) |> dirname |> dirname, "docs", "src", # hide
  "benchmark_python_bseries.txt")                                         # hide
try                                                                       # hide
  results = "```\n" * read(filename, String) * "```\n"                    # hide
  Markdown.parse(results)                                                 # hide
catch                                                                     # hide
  if get(ENV, "CI", nothing) != "true"                                    # hide
    println("We are not running CI so we do not show results here.")      # hide
  else                                                                    # hide
    rethrow()                                                             # hide
  end                                                                     # hide
end                                                                       # hide
````


Next, we look at the Python package
[`pybs`](https://github.com/henriksu/pybs)
and the following benchmark script. Note that this package does not provide
functionality for modifying integrators.

````@example
using BSeries, Markdown                                                   # hide
filename = joinpath(pathof(BSeries) |> dirname |> dirname, "docs", "src", # hide
  "benchmark_python_pybs.py")                                             # hide
script = "```python\n"                                                    # hide
for line in Iterators.drop(readlines(filename), 4)                        # hide
  startswith(line, "with") && continue                                    # hide
  line = replace(line, "  print" => "print")                              # hide
  global script                                                           # hide
  script = script * replace(line, ", file=io" => "") * "\n"               # hide
end                                                                       # hide
script = script * "```\n"                                                 # hide
Markdown.parse(script)                                                    # hide
````

The results are as follows.

````@example
using BSeries, Markdown                                                   # hide
filename = joinpath(pathof(BSeries) |> dirname |> dirname, "docs", "src", # hide
  "benchmark_python_pybs.txt")                                            # hide
try                                                                       # hide
  results = "```\n" * read(filename, String) * "```\n"                    # hide
  Markdown.parse(results)                                                 # hide
catch                                                                     # hide
  if get(ENV, "CI", nothing) != "true"                                    # hide
    println("We are not running CI so we do not show results here.")      # hide
  else                                                                    # hide
    rethrow()                                                             # hide
  end                                                                     # hide
end                                                                       # hide
````


Finally, we perform the same task using
[BSeries.jl](https://github.com/ranocha/BSeries.jl)
in Julia.

```@example
using BSeries, StaticArrays

A = @SArray [0 0; 1//2 0]
b = @SArray [0, 1//1]
c = @SArray [0, 1//2]
up_to_order = 9


println("Modified equation")
@time begin
  series = modified_equation(A, b, c, up_to_order)
  println(sum(values(series)))
end

@time begin
  series = modified_equation(A, b, c, up_to_order)
  println(sum(values(series)))
end


println("\nModifying integrator")
@time begin
  series = modifying_integrator(A, b, c, up_to_order)
  println(sum(values(series)))
end

@time begin
  series = modifying_integrator(A, b, c, up_to_order)
  println(sum(values(series)))
end
```



## References

[^RanochaKetcheson2020]:
  Hendrik Ranocha and David Ketcheson (2020)
  Energy Stability of Explicit Runge-Kutta Methods for Nonautonomous or
  Nonlinear Problems.
  SIAM Journal on Numerical Analysis
  [DOI: 10.1137/19M1290346](https://doi.org/10.1137/19M1290346)
