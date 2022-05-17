# [Symbolic computations](@id tutorial-symbolic-computations)

This tutorial describes some possibilities for symbolic computations based on
[BSeries.jl](https://github.com/ranocha/BSeries.jl). Currently, symbolic
computations in [BSeries.jl](https://github.com/ranocha/BSeries.jl) support
at least

- [SymEngine.jl](https://github.com/symengine/SymEngine.jl),
- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), and
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

as symbolic backends. You can find some performance comparisons of them in the
[benchmarks](@ref benchmarks-symbolic-computations).


## Generating LaTeX code

BSeries.jl allow fully symbolic computations of [`modified_equation`](@ref)s
and [`modifying_integrator`](@ref)s. More information about these features
are available in the respective tutorials. Here, we will also demonstrate how
LaTeX code can be generated via
[Latexify.jl](https://github.com/korsbo/Latexify.jl).

At the core, BSeries.jl is based on the representation of rooted trees via
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl). Thus, you need to
use the setup described in the docstring of `RootedTrees.latexify`.


### LaTeX code of a general modified equation

An explicit two-stage, second-order Runge-Kutta method is given by the Butcher
tableau

```math
\begin{array}{c|cc}
  0                 & 0                 &   \\
  \frac{1}{2\alpha} & \frac{1}{2\alpha} & 0 \\
  \hline
                    & 1 - \alpha        & \alpha\\
\end{array}
```

with real parameter ``\alpha``. We can compute the B-series of its modified
equation symbolically as follows.

```@example modified-equation-symengine
using BSeries, SymEngine

α = SymEngine.symbols("α")
A = [0 0; 1/(2α) 0]
b = [1-α, α]
c = [0, 1/(2α)]

series = modified_equation(A, b, c, 3)
```

To generate the appropriate LaTeX code, we need to remember that the B-series
coefficients shown above represent the modified vector field multiplied by the
time step size. Thus, we need to pass `reduce_order_by = 1` to `latexify`.
Other keyword arguments implemented for B-series are

- `f = "f"`: The symbol used for the indices of the elementary differentials.
- `dt = "h"`: The symbol used for the time step size. Alternatively, symbolic
  variables can also be used.
- `reduce_order_by = 0`: Reduce the power of the time step size accordingly.

```@example modified-equation-symengine
using Latexify
latexify(series, reduce_order_by=1, dt=SymEngine.symbols("h"), cdot=false) |> println
```


We can also use other packages for the symbolic computations, of course.
SymPy.jl often provides very clean expressions.

```@example modified-equation-sympy
using BSeries, SymPy

α = SymPy.symbols("α", real=true)
A = [0 0; 1/(2α) 0]
b = [1-α, α]
c = [0, 1/(2α)]

series = modified_equation(A, b, c, 3)
```

We can also generate LaTeX code as follows, using the same approach as for
SymEngine.jl.

```@example modified-equation-sympy
using Latexify
latexify(series, reduce_order_by=1, dt=SymPy.symbols("h"), cdot=false) |> println
```


Alternatively, we can also use Symbolics.jl.

```@example modified-equation-symbolics
using BSeries, Symbolics

Symbolics.@variables α h
A = [0 0; 1/(2α) 0]
b = [1-α, α]
c = [0, 1/(2α)]

series = modified_equation(A, b, c, 3)
```

## Working with rational coefficients

High-order Runge-Kutta methods are often given in terms of rational coefficients that have large denominators.
It's often best to use the rational coefficients in calculations, rather than converting to floating-point, since otherwise 
terms that ought to cancel may not cancel exactly.  Rational coefficients can be entered by using `//` for the division operation.

By default, Julia uses 64-bit integers for the numerator and denominator of a rational on 64-bit systems.  In practical 
calculations with high-order RK methods, the denominators may become too large to be represented with 64 bits, 
leading to overflow.  This can be remedied by specifying higher precision when entering the coefficients:

```@example int128-coefficients-symbolics
b = Rational{Int128}[(1081252805//134140608),(2639189439//74522560),(33646441//4191894),(-7873511875//210792384),(-504040617//14904512),(2110843561//115277085),(13//7),(1//2)];
```

If even `int128` is not enough, one can specify the type `BigInt`, which has adjustable precision.
