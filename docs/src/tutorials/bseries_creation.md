# Creating B-series

We have already seen some ways of creating B-series in the
[basic tutorial](@ref tutorial-bseries-basics). In this tutorial, we look at
ways to obtain the B-series of different time integration methods.


## B-series for Runge-Kutta methods

[BSeries.jl](https://github.com/ranocha/BSeries.jl) and
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl) provide
the type `RungeKuttaMethod` as wrapper of Butcher coefficients `A, b, c` of
Runge-Kutta methods.

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
to [`bseries`](@ref), you can also pass the wrapper struct. In fact, this is
the preferred method.

```@example ex:RK4
series = bseries(rk, 5)
```

We can check that the classical Runge-Kutta method is indeed fourth-order accurate.

```@example ex:RK4
series - ExactSolution(series)
```


## B-series for additive Runge-Kutta methods

[BSeries.jl](https://github.com/ranocha/BSeries.jl) and
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl) also support additive
Runge-Kutta methods via the wrapper `AdditiveRungeKuttaMethod`. For example,
we can write the Störmer-Verlet method as additive Runge-Kutta method following
Table II.2.1 of Hairer, Lubich, and Wanner (2002).

```@example ex:SV
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


### References

- Ernst Hairer, Gerhard Wanner, Christian Lubich.
  Geometric numerical integration.
  Springer, 2002.


## B-series for the average vector field method

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
\begin{align}
b(.) &= 1, \\
b([t_1, ..., t_n]) &= b(t_1)...b(t_n) / (n + 1).
\end{align}
```

by Celledoni, McLachlan, McLaren, Owren, Quispel, and Wright (2009). We can
implement this up to order 5 in BSeries.jl as follows.

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

We can check that this method is second-order accurate by comparing it to
the B-series of the exact solution, truncated at the same order.

```@example ex:AVF
series - ExactSolution(series)
```

### References

- Robert I. McLachlan, G. Reinout W. Quispel, and Nicolas Robidoux.
  "Geometric integration using discrete gradients."
  Philosophical Transactions of the Royal Society of London.
  Series A: Mathematical, Physical and Engineering Sciences 357,
  no. 1754 (1999): 1021-1045.
- G. Reinout W. Quispel, and David Ian McLaren.
  "A new class of energy-preserving numerical integration methods."
  Journal of Physics A: Mathematical and Theoretical 41, no. 4 (2008): 045206.
- Elena Celledoni, Robert I. McLachlan, David I. McLaren, Brynjulf Owren,
  G. Reinout W. Quispel, and William M. Wright.
  "Energy-preserving Runge-Kutta methods."
  ESAIM: Mathematical Modelling and Numerical Analysis 43, no. 4 (2009): 645-649.
  [DOI: 10.1051/m2an/2009020](https://doi.org/10.1051/m2an/2009020)
