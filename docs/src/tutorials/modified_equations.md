# [Modified equations](@id tutorial-modified-equation)

This tutorial describes the API of
[BSeries.jl](https://github.com/ranocha/BSeries.jl)
related to the notion of *modified equations*. The main API entry point is
[`modified_equation`](@ref).

Given a first-order autonomous ordinary differential equation (ODE)

```math
u'(t) = f(u(t))
```

and a B-series time integration method, the idea is to interpret the numerical
solution with given time step size ``h`` of the original ODE as the exact
solution of the modified ODE

```math
u'(t) = f_h(u(t)),
```

see for example [^ChartierHairerVilmart2010].


## Lotka-Volterra model

Here, we reproduce the example on p. 340 of [^HairerLubichWanner2006].
Thus, we consider the explicit Euler method to solve the classical
Lotka-Volterra model

```math
p'(t) = (2 - q) p,
\quad
q'(t) = (p - 1) q.
```

First, we set up the ODE and compute some numerical solutions using
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).

```@example ex:lotka-volterra
using OrdinaryDiffEq

function f(du, u, params, t)
  p, q = u
  dp = (2 - q) * p
  dq = (p - 1) * q
  du[1] = dp; du[2] = dq
  return nothing
end

u0 = [1.5, 2.25]
tspan = (0.0, 15.0)
ode = ODEProblem(f, u0, tspan)

dt = 0.1
sol_euler = solve(ode, Euler(), dt=dt)
sol_ref = solve(ode, Tsit5())
nothing # hide
```

Next, we look at some phase space plots of the numerical solution.

```@example ex:lotka-volterra
using LaTeXStrings, Plots

fig = plot(xguide=L"$q$", yguide=L"$p$")
default(linewidth=2)
plot!(fig, sol_ref, vars=(2, 1), label="Reference solution")
scatter!(fig, last.(sol_euler.u), first.(sol_euler.u),
         label="Explicit Euler, dt = $dt")
plot!(fig, xlims=(0.0, 9.0), ylims=(0.0, 5.0))

savefig(fig, "lotka_volterra_original.svg"); nothing # hide
```

![](lotka_volterra_original.svg)

The exact solution of this problem is periodic, but the explicit Euler method
produces an unstable trajectory. Here, we used an especially large time step to
more clearly illustrate what will follow, but the qualitative behavior is the
same for any time step size.

Next, we will derive the "modified equation" of the explicit Euler method and
solve this new ODE to high accuracy. The perturbed system takes the form of a
power series in the time step size `dt`, and in order to compute with it we will
truncate it at a certain order.

Here, we use [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)
for the symbolic computations.

```@example ex:lotka-volterra
using BSeries, StaticArrays, Symbolics

function solve_modified_equation(ode, truncation_orders, dt)
  # Explicit Euler method
  A = @SArray [0//1;]
  b = @SArray [1//1]
  c = @SArray [0//1]

  # Setup of symbolic variables
  @variables dt_sym
  u_sym = @variables p q
  f_sym = similar(u_sym); f(f_sym, u_sym, nothing, nothing)

  sol_euler = solve(ode, Euler(), dt=dt)

  fig = plot(xguide=L"$q$", yguide=L"$p$")
  default(linewidth=2)
  scatter!(fig, last.(sol_euler.u), first.(sol_euler.u),
          label="Explicit Euler, dt = $dt")

  for truncation_order in truncation_orders
    series = modified_equation(f_sym, u_sym, dt_sym, A, b, c, truncation_order)
    series = Symbolics.substitute.(series, dt_sym => dt)
    modified_f, _ = build_function(series, u_sym, expression=Val(false))
    modified_ode = ODEProblem((u, params, t) -> modified_f(u), ode.u0, tspan)
    modified_sol = solve(modified_ode, Tsit5())
    plot!(fig, modified_sol, vars=(2, 1),
          label="Modified ODE, order $(truncation_order-1)")
  end
  fig
end

fig = solve_modified_equation(ode, 2, dt)

savefig(fig, "lotka_volterra_modified1.svg"); nothing # hide
```

![](lotka_volterra_modified1.svg)

The exact solution of the Lotka-Volterra model is periodic, but Euler's method
generates a solution with growing amplitude. The modified equations accurately
predict this.

Now we go to the next order and increase the time step size `dt` slightly.

```@example ex:lotka-volterra
fig = solve_modified_equation(ode, 2:3, 0.11)

savefig(fig, "lotka_volterra_modified2.svg"); nothing # hide
```

![](lotka_volterra_modified2.svg)

Using a larger step size, we see that the first-order modified equations are
not fully accurate, but by including the ``O(h^2)`` terms we get much better
accuracy at late times. Let's keep going.

```@example ex:lotka-volterra
fig = solve_modified_equation(ode, 2:4, 0.12)

savefig(fig, "lotka_volterra_modified3.svg"); nothing # hide
```

![](lotka_volterra_modified3.svg)


## References

[^HairerLubichWanner2006]:
  Ernst Hairer, Christian Lubich, Gerhard Wanner (2006)
  Geometric Numerical Integration.
  [DOI: 10.1007/3-540-30666-8](https://doi.org/10.1007/3-540-30666-8)

[^ChartierHairerVilmart2010]:
  Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
