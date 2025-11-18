# [Creating B-series](@id tutorial-bseries-creation)

We have already seen some ways of creating B-series in the
[basic tutorial](@ref tutorial-bseries-basics). In this tutorial, we look at
ways to obtain the B-series of different time integration methods.


## B-series for Runge-Kutta methods

[BSeries.jl](https://github.com/ranocha/BSeries.jl) and
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl) provide
the type [`RungeKuttaMethod`](@ref BSeries.RootedTrees.RungeKuttaMethod)
as wrapper of Butcher coefficients `A, b, c` of Runge-Kutta methods.
For example, you can create the classical explicit, fourth-order Runge-Kutta
method as follows.

```@example ex:RK4
using BSeries

# classical RK4
A = [0 0 0 0
      1//2 0 0 0
      0 1//2 0 0
      0 0 1 0]
b = [1 // 6, 1 // 3, 1 // 3, 1 // 6]
rk = RungeKuttaMethod(A, b)
```

Instead of passing the Butcher coefficients explicitly
to [`bseries`](@ref), you can also pass the wrapper struct `rk`. In fact,
this is the preferred method.

```@example ex:RK4
series = bseries(rk, 5)
```

We can check that the classical Runge-Kutta method is indeed fourth-order accurate.

```@example ex:RK4
series - ExactSolution(series)
```

```@example ex:RK4
order_of_accuracy(series)
```


## [B-series for Runge-Kutta methods via composition](@id tutorial-bseries-creation-RK-compose)

The method described in the previous section is the preferred way of constructing
B-series for Runge-Kutta methods. However, it is also possible to construct
the B-series of a Runge-Kutta method manually via composition of simpler B-series.
First, we start with the B-series of the identity interpreted as a time integrator
``u^{n+1} = u^n``. This B-series is called the [`IdentityMap`](@ref) in
[BSeries.jl](https://github.com/ranocha/BSeries.jl).

```@example ex:RK2-compose
using BSeries
y1 = bseries(IdentityMap{Int}(), 4)
```

This will be the first stage of a two-stage explicit Runge-Kutta method with
symbolic coefficients. The next stage is given by
``y^2 = y^1 + a_{21} \Delta t f(y^1)``.

First, we create the B-series of ``\Delta t f(y^1)`` by composing
the B-series of the first stage with the B-series of the vector field ``h f``,
see [`compose`](@ref) and [`IdentityField`](@ref).

```@example ex:RK2-compose
compose(y1, IdentityField())
```

Then, we can compute the second stage.

```@example ex:RK2-compose
using SymPyPythonCall
a21 = symbols("a21", real = true)
y2 = y1 + a21 * compose(y1, IdentityField())
```

Similarly, we can compute the new update
``u^{n+1} = u^n + b_1 \Delta t f(y^1) + b_2 \Delta t f(y^2)``.

```@example ex:RK2-compose
b1, b2 = symbols("b1 b2", real = true)
series = y1 + b1 * compose(y1, IdentityField()) + b2 * compose(y2, IdentityField())
```

Alternatively, we can obtain the same result by using the
[`RungeKuttaMethod`](@ref) wrapper as shown in the previous section.

```@example ex:RK2-compose
A = [0 0; a21 0]
b = [b1, b2]
rk = RungeKuttaMethod(A, b)
series_rk = bseries(rk, order(series))
@assert series == series_rk # hide
```


## B-series for additive Runge-Kutta methods

[BSeries.jl](https://github.com/ranocha/BSeries.jl) and
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl) also support additive
Runge-Kutta methods via the wrapper
[`AdditiveRungeKuttaMethod`](@ref BSeries.RootedTrees.AdditiveRungeKuttaMethod).
For example, we can write the Störmer-Verlet method as an additive Runge-Kutta
method following Table II.2.1 of Hairer, Lubich, and Wanner (2002).

```@example ex:SV
using BSeries

As = [
    [0 0; 1//2 1//2],
    [1//2 0; 1//2 0],
]
bs = [
    [1 // 2, 1 // 2],
    [1 // 2, 1 // 2],
]
ark = AdditiveRungeKuttaMethod(As, bs)
```

We can create the B-series as usual, truncated to order 3.

```@example ex:SV
series = bseries(ark, 3)
```

For additive Runge-Kutta methods like this, we use colored rooted trees.
Again, we can check the order of accuracy by comparing the coefficients to
the exact solution.

```@example ex:SV
series - ExactSolution(series)
```

```@example ex:SV
order_of_accuracy(series)
```

We can also create LaTeX code for this B-series as follows.

```@example ex:SV
using Latexify
latexify(series)
```

When compiled using the preamble code shown in the docstring of
`RootedTrees.latexify`, the output looks as follows.

![bseries-SV](https://user-images.githubusercontent.com/12693098/216765237-0ec1138f-e3bb-45e4-bc81-9de6a3904e46.jpg)


### References

- Ernst Hairer, Gerhard Wanner, Christian Lubich.
  Geometric numerical integration.
  Springer, 2002.
  [DOI: 10.1007/3-540-30666-8](https://doi.org/10.1007/3-540-30666-8)


## B-series for Rosenbrock methods

[BSeries.jl](https://github.com/ranocha/BSeries.jl) and
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl) also support
Rosenbrock (Rosenbrock-Wanner, ROW) methods via the wrapper
[`RosenbrockMethod`](@ref BSeries.RootedTrees.RosenbrockMethod).
For example, a classical ROW method of Kaps and Rentrop (1979) can be
parameterized as follows.

```@example ex:ROW
using BSeries

γ = [0.395 0 0 0;
     -0.767672395484 0.395 0 0;
     -0.851675323742 0.522967289188 0.395 0;
     0.288463109545 0.880214273381e-1 -0.337389840627 0.395]
A = [0 0 0 0;
     0.438 0 0 0;
     0.796920457938 0.730795420615e-1 0 0;
     0.796920457938 0.730795420615e-1 0 0]
b = [0.199293275701, 0.482645235674, 0.680614886256e-1, 0.25]
ros = RosenbrockMethod(γ, A, b)
```

We can create the B-series as usual, truncated to order 5.

```@example ex:ROW
series = bseries(ros, 5)
```

Again, we can check the order of accuracy by comparing the coefficients to
the exact solution.

```@example ex:ROW
series - ExactSolution(series)
```

```@example ex:ROW
order_of_accuracy(series)
```

### References

- Peter Kaps and Peter Rentrop.
  "Generalized Runge-Kutta methods of order four with stepsize control for
  stiff ordinary differential equations."
  Numerische Mathematik 33, no. 1 (1979): 55-68.
  [DOI: 10.1007/BF01396495](https://doi.org/10.1007/BF01396495)


## [B-series for the average vector field method](@id tutorial-bseries-creation-AVF)

Consider the autonomous ODE

```math
u'(t) = f\bigl( u(t) \bigr).
```

The average vector field method

```math
u^{n+1} = u^{n} + \Delta t \int_0^1 f\bigl(\xi u^{n+1} + (1 - \xi) u^{n}\bigr) \mathrm{d} \xi
```

introduced by McLachlan, Quispel, and Robidoux (1999) is a second-order accurate
method. Quispel and McLaren (2008) discovered that it is indeed a B-series method.
Its coefficients are given explicitly as

```math
\begin{align*}
b(.) &= 1, \\
b([t_1, ..., t_n]) &= b(t_1)...b(t_n) / (n + 1).
\end{align*}
```

by Celledoni, McLachlan, McLaren, Owren, Quispel, and Wright (2009). We can
implement this up to order 5 in [BSeries.jl](https://github.com/ranocha/BSeries.jl)
as follows.

```@example ex:AVF
using BSeries

series = bseries(5) do t, series
    if order(t) in (0, 1)
        return 1 // 1
    else
        v = 1 // 1
        n = 0
        for subtree in SubtreeIterator(t)
            v *= series[subtree]
            n += 1
        end
        return v / (n + 1)
    end
end
```

[BSeries.jl](https://github.com/ranocha/BSeries.jl) also offers a
convenience constructor using the type [`AverageVectorFieldMethod`](@ref)
as follows.

```@example ex:AVF
series == bseries(AverageVectorFieldMethod(), 5)
```

We can check that this method is second-order accurate by comparing it to
the B-series of the exact solution, truncated at the same order.

```@example ex:AVF
series - ExactSolution(series)
```

```@example ex:AVF
order_of_accuracy(series)
```

### References

- Robert I. McLachlan, G. Reinout W. Quispel, and Nicolas Robidoux.
  "Geometric integration using discrete gradients."
  Philosophical Transactions of the Royal Society of London.
  Series A: Mathematical, Physical and Engineering Sciences 357,
  no. 1754 (1999): 1021-1045.
  [DOI: 10.1098/rsta.1999.0363](https://doi.org/10.1098/rsta.1999.0363)
- G. Reinout W. Quispel, and David Ian McLaren.
  "A new class of energy-preserving numerical integration methods."
  Journal of Physics A: Mathematical and Theoretical 41, no. 4 (2008): 045206.
  [DOI: 10.1088/1751-8113/41/4/045206](https://doi.org/10.1088/1751-8113/41/4/045206)
- Elena Celledoni, Robert I. McLachlan, David I. McLaren, Brynjulf Owren,
  G. Reinout W. Quispel, and William M. Wright.
  "Energy-preserving Runge-Kutta methods."
  ESAIM: Mathematical Modelling and Numerical Analysis 43, no. 4 (2009): 645-649.
  [DOI: 10.1051/m2an/2009020](https://doi.org/10.1051/m2an/2009020)
