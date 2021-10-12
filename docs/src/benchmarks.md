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
Pkg.status(["SymEngine", "SymPy", "Symbolics"])
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
simple performance comparisons. Again, we just use `@time` twice to get an idea
of the performance after compilation, allowing us to compare orders of magnitude.

First, we start with the Python package
[`BSeries`](https://github.com/ketch/BSeries).
```@example benchmark-Python-BSeries
using PyCall

py"""
import BSeries.bs as bs
import nodepy.runge_kutta_method as rk

midpoint_method = rk.loadRKM("Mid22")
up_to_order = 9
"""

@time py"""
series = bs.modified_equation(None, None,
                              midpoint_method.A, midpoint_method.b,
                              up_to_order, True)
print(sum(series.values()))
"""

@time py"""
series = bs.modified_equation(None, None,
                              midpoint_method.A, midpoint_method.b,
                              up_to_order, True)
print(sum(series.values()))
"""
```


Next, we look at the Python package [`pybs`](https://github.com/henriksu/pybs).
```@example
using PyCall
println(PyCall.python)

py"""
import sys
print(sys.version_info)
"""

PyCall.Conda.list()
```


```@example
using PyCall
PyCall.python_cmd(`-c "import pybs; print(pybs.__path__)"`) |> run
```


```@example benchmark-Python-pybs
using PyCall

py"""
import pybs
"""
```

```@example benchmark-Python-pybs
using PyCall

py"""
import pybs
print(pybs.__path__)
"""
```

```@example benchmark-Python-pybs
using PyCall

py"""
import pybs
from pybs.rungekutta import methods as rk_methods

midpoint_method = rk_methods.RKmidpoint
up_to_order = 9
number_of_terms = pybs.unordered_tree.number_of_trees_up_to_order(up_to_order+1)

from itertools import islice
def first_values(f, n):
  return (f(tree) for tree in islice(pybs.unordered_tree.tree_generator(), 0, n))

"""
```

```@example benchmark-Python-pybs
@time py"""
midpoint_series = midpoint_method.phi()
# series = pybs.series.modified_equation(midpoint_series)
# print(sum(first_values(series, number_of_terms)))
"""

```@example benchmark-Python-pybs
@time py"""
midpoint_series = midpoint_method.phi()
series = pybs.series.modified_equation(midpoint_series)
# print(sum(first_values(series, number_of_terms)))
"""

```@example benchmark-Python-pybs
@time py"""
midpoint_series = midpoint_method.phi()
series = pybs.series.modified_equation(midpoint_series)
print(sum(first_values(series, number_of_terms)))
"""

@time py"""
midpoint_series = midpoint_method.phi()
series = pybs.series.modified_equation(midpoint_series)
print(sum(first_values(series, number_of_terms)))
"""
```


Finally, we perform the same task using
[BSeries.jl](https://github.com/ranocha/BSeries.jl).
```@example
using BSeries, StaticArrays

A = @SArray [0 0; 1//2 0]
b = @SArray [0, 1//1]
c = @SArray [0, 1//2]
up_to_order = 9

@time begin
  series = modified_equation(A, b, c, up_to_order)
  println(sum(values(series)))
end

@time begin
  series = modified_equation(A, b, c, up_to_order)
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
