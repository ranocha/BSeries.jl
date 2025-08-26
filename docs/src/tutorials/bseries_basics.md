# [B-series basics](@id tutorial-bseries-basics)

In this tutorial we use BSeries.jl to investigate error expansions for Runge-Kutta (RK)
methods, generically or when applied to a specific ordinary differential equation (ODE).


```@example bseries-basics
# Load the packages we will use.
# These must first be installed using: import Pkg; Pkg.add("package_name")
using BSeries
using Latexify  # Only needed for some pretty-printing cells below using `latexify`
import SymPyPythonCall; sp = SymPyPythonCall;
```

## B-series for a generic ODE

First we specify the Butcher coefficients of the RK method.
This can include symbolic expressions and parameterized families of methods.
Here is a generic 2-stage, 2nd-order method:


```@example bseries-basics
α = sp.symbols("α", real = true)
A = [0 0; 1/(2*α) 0]; b = [1-α, α]; c = [0, 1/(2*α)]
coeffs2 = bseries(A, b, c, 3)
latexify(coeffs2, cdot=false)
```

The rooted trees are printed as nested lists, essentially in the form used in
Butcher's book. The rooted trees written in this way can be rendered in LaTeX
using the package [`forest`](https://ctan.org/pkg/forest); unfortunately,
there is no easy way to render them in the browser. Nevertheless, you can render
them using LaTeX with an appropriate preamble, see the docstring of
`RootedTrees.latexify`. The rendered output looks like this:

![bseries_creation_eq1-1](https://user-images.githubusercontent.com/12693098/193994163-e53d24a8-f74e-4f95-b07d-225ebde83f70.png)

We have generated the B-series up to terms of order $h^3$.  The terms $F_f()$
represent elementary differentials, which are products of derivatives of the
ODE right-hand side.  Since we haven't specified an ODE, these are indicated
simply by the associated rooted tree.


Here is a B-series for the classical 4th-order method, expanded up to 5th-order terms:

```@example bseries-basics
A = [0 0 0 0; 1//2 0 0 0; 0 1//2 0 0; 0 0 1 0];
b = [1//6, 1//3, 1//3, 1//6];
c = [0, 1//2, 1//2, 1];

coeffs4 = bseries(A, b, c, 5)
latexify(coeffs4, cdot=false)
```

![bseries_creation_eq2-1](https://user-images.githubusercontent.com/12693098/193994166-a9178001-702d-4f9b-a3f6-6a89251ddb7f.png)

We can also print out the B-series coefficients this way:


```@example bseries-basics
coeffs4
```

In this form, the rooted trees are printed as level sequences.
The corresponding coefficients are on the right.

You can use the function `RootedTrees.set_printing_style` to change the
printing style globally. For example, you can use the notation of Butcher
as follows.

```@example bseries-basics
RootedTrees.set_printing_style("butcher")
coeffs4
```

To use the level sequence representation, you need to change the printing style
again.

```@example bseries-basics
RootedTrees.set_printing_style("sequence")
coeffs4
```


## Exact series and local error

We can also get the B-series of the exact solution:


```@example bseries-basics
coeffs_ex = ExactSolution(coeffs4)
```


```@example bseries-basics
latexify(coeffs_ex, cdot=false)
```

![bseries_creation_eq3-1](https://user-images.githubusercontent.com/12693098/193994175-22356d01-edb9-44b6-afd3-4354a3daffc6.png)

We can find the local error by subtracting the exact solution B-series from the RK method B-series:

```@example bseries-basics
latexify(coeffs4-coeffs_ex, cdot=false)
```

![bseries_creation_eq4-1](https://user-images.githubusercontent.com/12693098/193994179-7ffcced2-6760-46fc-829a-d6c5814d543f.png)

This confirms that the method is of 4th order, since all terms involving
smaller powers of $h$ vanish exactly.  We don't see the $h^6$ and higher
order terms since we only generated the truncated B-series up to 5th order.
We can also obtain the order of accuracy without comparing the coefficients
to the exact solution manually:

```@example bseries-basics
order_of_accuracy(coeffs4)
```

For the 2nd-order method, we get

```@example bseries-basics
order_of_accuracy(coeffs2)
```

with the following leading error terms:

```@example bseries-basics
latexify(coeffs2-coeffs_ex, cdot=false)
```

![bseries_creation_eq5-1](https://user-images.githubusercontent.com/12693098/193994181-108aa3a7-e2fb-4247-a770-9647ebe861c8.png)

This confirms again the accuracy of the method, and shows us that we
can eliminate one of the leading error terms completely if we take
$\alpha=3/4$ (this is known as Ralston's method, or sometimes as Heun's method).


## B-series for a specific ODE

Next, let us define an ODE.  We'll consider the Prothero-Robinson problem:

```math
    y'(t) = \lambda(y-\sin(t)) + \cos(t).
```

For a non-autonomous ODE like this, it's convenient to rewrite the problem
in autonomous form.  We set $u=[y,t]^T$ and

```math
\begin{align}
u_1'(t) & = \lambda(u_1 - \sin(u_2)) + \cos(u_2) \\
u_2'(t) & = 1.
\end{align}
```


```@example bseries-basics
λ = sp.symbols("λ", real = true)
y, t = sp.symbols("y t", real = true)
h = sp.symbols("h", real = true)

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
