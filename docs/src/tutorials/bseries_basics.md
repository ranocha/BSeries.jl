# [B-Series basics](@id tutorial-bseries-basics)
In this tutorial we use `bseries.jl` to investigate error expansions for RK
methods, generically or when applied to a specific ODE.


```@example bseries-basics
# Load the packages we will use.  These must first be installed using: import Pkg; Pkg.add("package_name")
using BSeries
using Latexify  # Only needed for some pretty-printing cells below
import SymPy; sp=SymPy;
```

# B-Series for a generic ODE

First we specify the Butcher coefficients of the RK method.
This can include symbolic expressions and parameterized families of methods.
Here is a generic 2-stage, 2nd-order method:


```@example bseries-basics
α = sp.symbols("α", real=true)
A = [0 0; 1/(2*α) 0]; b = [1-α, α]; c = [0, 1/(2*α)]
coeffs2 = bseries(A,b,c,3)
latexify(coeffs2, cdot=false)
```

We have generated the B-series up to terms of order $h^3$.  The terms $F_f()$
represent elementary differentials, which are products of derivatives of the
ODE right-hand side.  Since we haven't specified an ODE, these are indicated
simply by the associated rooted tree.  The rooted trees are printed as nested
lists, essentially in the form used in Butcher's book.  The rooted trees written
in this way can be rendered in LaTex using the package `forest`; unfortunately,
there is no easy way to render them in the browser.

Here is a B-series for a 4th-order method, expanded up to 5th-order terms:


```@example bseries-basics
A = Rational{Int128}[0 0 0 0 0 0 0 0;(-1//6) (1//2) 0 0 0 0 0 0;(-1//10) (1//10) (1//2) 0 0 0 0 0;(-21463//39375) (21017//26250) (-5//9) (1//2) 0 0 0 0;(-59588//54675) (118717//36450) (-4375//2187) 0 (1//2) 0 0 0;(-19993033//9443328) (28508695//3147776) (-13577105//2360832) (-4090625//3147776) (1136025//3147776) (1//2) 0 0;(367020141781//199294617600) (814214904871//22143846400) (-29834937659//1992946176) (-1983358776875//87689631744) (-6702625935//885753856) (688576//109395) (1//2) 0;(1081252805//134140608) (2639189439//74522560) (33646441//4191894) (-7873511875//210792384) (-504040617//14904512) (2110843561//115277085) (13//7) (1//2)];
b = Rational{Int128}[(1081252805//134140608),(2639189439//74522560),(33646441//4191894),(-7873511875//210792384),(-504040617//14904512),(2110843561//115277085),(13//7),(1//2)];
c = Rational{Int128}[0,(1//3),(1//2),(1//5),(2//3),(3//4),(1//4),1];

coeffs4 = bseries(A,b,c,5)
latexify(coeffs4, cdot=false)
```

We can also print out the B-series coefficients this way:


```@example bseries-basics
coeffs4
```

In this form, the rooted trees are printed as level sequences.  The corresponding coefficients are on the right.

# Exact series and local error

We can also get the B-series of the exact solution:


```@example bseries-basics
coeffs_ex = ExactSolution(coeffs4)
```


```@example bseries-basics
latexify(coeffs_ex,cdot=false)
```

We can find the local error by subtracting the exact solution B-series from the RK method B-series:


```@example bseries-basics
latexify(coeffs4-coeffs_ex,cdot=false)
```


This confirms that the method is of 4th order, since all terms involving
smaller powers of $h$ vanish exactly.  We don't see the $h^6$ and higher
order terms since we only generated the truncated B-series up to 5th order.

For the 2nd-order method, we get:


```@example bseries-basics
latexify(coeffs2-coeffs_ex,cdot=false)
```

This confirms again the accuracy of the method, and shows us that we
can eliminate one of the leading error terms completely if we take
$\alpha=3/4$ (this is known as Ralston's method, or sometimes as Heun's method).

# B-Series for a specific ODE

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
evaluate(ff,u,h,coeffs4)[1]
```

Notice that the series is truncated at the same order that we specified
when we initially generated it from the RK coefficients.

Here's the B-Series for the exact solution of the same ODE:


```@example bseries-basics
evaluate(ff,u,h,coeffs_ex)[1]
```

And their difference, which is the local error:


```@example bseries-basics
expr = sp.simplify(evaluate(ff,u,h,coeffs4)-evaluate(ff,u,h,coeffs_ex))[1]
```

Because of the simplicity of the PR problem, the error has a relatively simple form:


```@example bseries-basics
sp.simplify(sp.expand(expr))
```

```@example bseries-basics
sp.collect(sp.expand(expr),λ)
```

# B-series for a generic RK method

We can also examine just the elementary differentials, without specifying a RK method:


```@example bseries-basics
elementary_differentials(ff,u,5)
```
