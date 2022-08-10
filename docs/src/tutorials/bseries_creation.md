# [Creating B-series](@id tutorial-bseries-creation)

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


### References

- Ernst Hairer, Gerhard Wanner, Christian Lubich.
  Geometric numerical integration.
  Springer, 2002.


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


## [Composing B-series](@id tutorial-bseries-creation-compose)

B-series representing a mapping such as a time integration method can be
[`compose`](@ref)d. This operation is equivalent to performing a step with
one method after another. Here, we present the example of Butcher's method of
effective order 5. This is a fourth-order Runge-Kutta method resulting in a
fifth-order method when composed with a special starting and finishing
procedure.

First, we set up the B-series of the main method (method "a" in Butcher's paper):

```@example ex:compose
using BSeries

A = [0 0 0 0 0;
     1//5 0 0 0 0;
     0 2//5 0 0 0;
     3//16 0 5//16 0 0;
     1//4 0 -5//4 2 0]
b = [1 // 6, 0, 0, 2 // 3, 1 // 6]
rk_a = RungeKuttaMethod(A, b)
series_a = bseries(rk_a, 6)
```

Note that this method is fourth-order accurate:

```@example ex:compose
series_a - ExactSolution(series_a)
```

Next, we set up the starting procedure (method "b" in Butcher's paper):

```@example ex:compose
A = [0 0 0 0 0;
     1//5 0 0 0 0;
     0 2//5 0 0 0;
     75//64 -9//4 117//64 0 0;
     -37//36 7//3 -3//4 4//9 0]
b = [19 // 144, 0, 25 // 48, 2 // 9, 1 // 8]
rk_b = RungeKuttaMethod(A, b)
series_b = bseries(rk_b, 6)
```

Note that this method is only third-order accurate - as is the finishing
procedure given by

```@example ex:compose
A = [0 0 0 0 0;
     1//5 0 0 0 0;
     0 2//5 0 0 0;
     161//192 -19//12 287//192 0 0;
     -27//28 19//7 -291//196 36//49 0]
b = [7 // 48, 0, 475 // 1008, 2 // 7, 7 // 72]
rk_c = RungeKuttaMethod(A, b)
series_c = bseries(rk_c, 6)
```

Finally, we can compose the three methods to obtain

```@example ex:compose
series = compose(series_b, series_a, series_c, normalize_stepsize = true)
```

Note that this composition has to be read from left to right. Finally, we check
that the resulting `series` is indeed fifth-order accurate:

```@example ex:compose
series - ExactSolution(series)
```


### References

- Butcher, J. C.
  "The effective order of Runge-Kutta methods."
  In Conference on the numerical solution of differential equations,
  pp. 133-139. Springer, Berlin, Heidelberg, 1969.
  https://doi.org/10.1007/BFb0060019
