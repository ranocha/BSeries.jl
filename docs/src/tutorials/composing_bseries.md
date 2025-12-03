# Composing B-series

B-series can be [`compose`](@ref)d, i.e., the result of one B-series
can be inserted into another B-series. We saw this operation already
in a [previous tutorial](@ref tutorial-bseries-creation-RK-compose)
creating the B-series of a Runge-Kutta method step by step. Here,
we present additional examples of composing B-series.


## [Richardson extrapolation](@id compose-richardson-extrapolation)

Richardson extrapolation is a well-known technique to increase the
order of accuracy of a time integration method. Assume we are given
a Runge-Kutta method with step size ``h``. Let ``u^{n+1}_h`` be the
solution after one step with step size ``h``, and let
``u^{n+1}_{2 \times \frac{h}{2}}`` be the solution after two steps with
step size ``\frac{h}{2}``. The difference between these two approximations
can be used to estimate the leading-order error term and adapt the
step size accordingly.

Moreover, we can combine the two approximations to obtain a new approximation
``u^{n+1}_R`` by extrapolating to step size ``h = 0`` as follows:

```math
u^{n+1}_R = u^{n+1}_{2 \times \frac{h}{2}} + \frac{ u^{n+1}_{2 \times \frac{h}{2}} - u^{n+1}_h }{ 2^p - 1 },
```

where ``p`` is the order of accuracy of the original method. First,
we use the classical fourth-order Runge-Kutta method as our base method.

```@example ex:richardson
using BSeries

# classical RK4
A = [0 0 0 0
      1//2 0 0 0
      0 1//2 0 0
      0 0 1 0]
b = [1 // 6, 1 // 3, 1 // 3, 1 // 6]
rk = RungeKuttaMethod(A, b)
series_rk = bseries(rk, 6)
```

```@example ex:richardson
order_of_accuracy(series_rk)
```

Next, we compute the B-series of two steps with half the step size.

```@example ex:richardson
series_rk_2 = compose(series_rk, series_rk, normalize_stepsize = true)
```

```@example ex:richardson
order_of_accuracy(series_rk_2)
```

Finally, we can combine the two B-series to obtain the Richardson-extrapolated
B-series.

```@example ex:richardson
p = order_of_accuracy(series_rk)
series_richardson = series_rk_2 + (series_rk_2 - series_rk) / (2^p - 1)
```

```@example ex:richardson
order_of_accuracy(series_richardson)
```

We can also perform this composition of the methods at the level of the
Butcher tableau.

```@example ex:richardson
AA = [A zero(A) zero(A)
      zero(A) A/2 zero(A)
      zero(A) repeat(b'/2, length(b)) A/2]
bb = vcat(-b / (2^p - 1),
          b / 2 * (1 + 1 // (2^p - 1)),
          b / 2 * (1 + 1 // (2^p - 1)))
rk_richardson = RungeKuttaMethod(AA, bb)
series_richardson_2 = bseries(rk_richardson, 6)
all(iszero, values(series_richardson - series_richardson_2))
```


## [Effective order and composition](@id compose-effective-order)

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
order_of_accuracy(series_a)
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

```@example ex:compose
order_of_accuracy(series_b)
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

```@example ex:compose
order_of_accuracy(series_c)
```

Finally, we can compose the three methods to obtain

```@example ex:compose
series = compose(series_b, series_a, series_c, normalize_stepsize = true)
```

Note that this composition has to be read from left to right. Finally, we check
that the resulting `series` is indeed fifth-order accurate:

```@example ex:compose
order_of_accuracy(series)
```

```@example ex:compose
series - ExactSolution(series)
```

### References

- Butcher, J. C.
  "The effective order of Runge-Kutta methods."
  In Conference on the numerical solution of differential equations,
  pp. 133-139. Springer, Berlin, Heidelberg, 1969.
  [DOI: 10.1007/BFb0060019](https://doi.org/10.1007/BFb0060019)

