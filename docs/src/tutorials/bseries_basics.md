# [B-Series basics](@id tutorial-bseries-basics)

In this tutorial we use BSeries.jl to investigate error expansions for Runge-Kutta (RK)
methods, generically or when applied to a specific ordinary differential equation (ODE).


```@example bseries-basics
# Load the packages we will use.
# These must first be installed using: import Pkg; Pkg.add("package_name")
using BSeries
using Latexify  # Only needed for some pretty-printing cells below using `latexify`
import SymPy; sp=SymPy;
```

## B-Series for a generic ODE

First we specify the Butcher coefficients of the RK method.
This can include symbolic expressions and parameterized families of methods.
Here is a generic 2-stage, 2nd-order method:


```@example bseries-basics
α = sp.symbols("α", real=true)
A = [0 0; 1/(2*α) 0]; b = [1-α, α]; c = [0, 1/(2*α)]
coeffs2 = bseries(A, b, c, 3)
latexify(coeffs2, cdot=false)
```

We have generated the B-series up to terms of order $h^3$.  The terms $F_f()$
represent elementary differentials, which are products of derivatives of the
ODE right-hand side.  Since we haven't specified an ODE, these are indicated
simply by the associated rooted tree.  The rooted trees are printed as nested
lists, essentially in the form used in Butcher's book.  The rooted trees written
in this way can be rendered in LaTeX using the package [`forest`](https://ctan.org/pkg/forest); unfortunately,
there is no easy way to render them in the browser.

Here is a B-series for the classical 4th-order method, expanded up to 5th-order terms:


```@example bseries-basics
A = [0 0 0 0; 1//2 0 0 0; 0 1//2 0 0; 0 0 1 0];
b = [1//6, 1//3, 1//3, 1//6];
c = [0, 1//2, 1//2, 1];

coeffs4 = bseries(A, b, c, 5)
latexify(coeffs4, cdot=false)
```

We can also print out the B-series coefficients this way:


```@example bseries-basics
coeffs4
```

In this form, the rooted trees are printed as level sequences.
The corresponding coefficients are on the right.


## Exact series and local error

We can also get the B-series of the exact solution:


```@example bseries-basics
coeffs_ex = ExactSolution(coeffs4)
```


```@example bseries-basics
latexify(coeffs_ex, cdot=false)
```

We can find the local error by subtracting the exact solution B-series from the RK method B-series:


```@example bseries-basics
latexify(coeffs4-coeffs_ex, cdot=false)
```


This confirms that the method is of 4th order, since all terms involving
smaller powers of $h$ vanish exactly.  We don't see the $h^6$ and higher
order terms since we only generated the truncated B-series up to 5th order.

For the 2nd-order method, we get:


```@example bseries-basics
latexify(coeffs2-coeffs_ex, cdot=false)
```

This confirms again the accuracy of the method, and shows us that we
can eliminate one of the leading error terms completely if we take
$\alpha=3/4$ (this is known as Ralston's method, or sometimes as Heun's method).


## B-Series for a specific ODE

Next, let us define an ODE.  We'll consider the Prothero-Robinson problem:

```math
    y'(t) = \lambda(y-\sin(t)) + \cos(t).
```

For a non-autonomous ODE like this, it's convenient to rewrite the problem
in autonomous form.  We set $u=[y,t]^T$ and

```math
\begin{align}
u_1'(t) & = \lambda(u_1 - \sin(u_2)) + \cos(u_2) \\
u_2'(t) & = t.
\end{align}
```


```@example bseries-basics
λ = sp.symbols("λ", real=true)
y, t = sp.symbols("y t", real=true)
h = sp.symbols("h", real=true)

u = [y, t]
ff = [λ*(u[1]-sin(t))+cos(t), 1]
```

Finally, we get the B-Series for our RK method applied to our ODE:


```@example bseries-basics
evaluate(ff, u, h, coeffs4)[1]
```

Notice that the series is truncated at the same order that we specified
when we initially generated it from the RK coefficients.

Here's the B-Series for the exact solution of the same ODE:


```@example bseries-basics
evaluate(ff, u, h, coeffs_ex)[1]
```

And their difference, which is the local error:


```@example bseries-basics
expr = sp.simplify(evaluate(ff, u, h, coeffs4) - evaluate(ff, u, h,coeffs_ex))[1]
```


## B-series for a generic RK method

We can also examine just the elementary differentials, without specifying a RK method:


```@example bseries-basics
elementary_differentials(ff, u, 5)
```


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
It's coefficients are given explicitly as

```math
\begin{align}
b(.) &= 1, \\
b([t_1, ..., t_n]) &= b(t_1)...b(t_n) / (n + 1).
```

by Celledoni. McLachlan, McLaren, Owren, Quispel, and Wright (2009). We can
implement this up to order 5 in BSeries.jl as follows.

```@example ex:AVF
using BSeries
series = bseries(5) do t, series
    if order(t) in (0, 1)
        return 1//1
    else
        v = 1//1
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
  "Energy-preserving runge-kutta methods."
  ESAIM: Mathematical Modelling and Numerical Analysis 43, no. 4 (2009): 645-649.
  [DOI: 10.1051/m2an/2009020](https://doi.org/10.1051/m2an/2009020)
