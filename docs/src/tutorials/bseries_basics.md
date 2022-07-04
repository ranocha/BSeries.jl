# [B-Series basics](@id tutorial-bseries-basics)
In this tutorial we use `bseries.jl` to investigate error expansions for RK
methods, generically or when applied to a specific ODE.


```julia
# Load the packages we will use.  These must first be installed using: import Pkg; Pkg.add("package_name")
using BSeries
using Latexify  # Only needed for some pretty-printing cells below
import SymPy; sp=SymPy;
```

# B-Series for a generic ODE

First we specify the Butcher coefficients of the RK method.
This can include symbolic expressions and parameterized families of methods.
Here is a generic 2-stage, 2nd-order method:


```julia
α = sp.symbols("α", real=true)
A = [0 0; 1/(2*α) 0]; b = [1-α, α]; c = [0, 1/(2*α)]
coeffs2 = bseries(A,b,c,3)
latexify(coeffs2, cdot=false)
```




$F_{f}\mathopen{}\left( \varnothing \right)\mathclose{} + h F_{f}\mathopen{}\left( \rootedtree[] \right)\mathclose{} + \frac{h^{2}}{2} F_{f}\mathopen{}\left( \rootedtree[[]] \right)\mathclose{} + \frac{h^{3}}{8 \alpha} F_{f}\mathopen{}\left( \rootedtree[[][]] \right)\mathclose{}$



We have generated the B-series up to terms of order $h^3$.  The terms $F_f()$
represent elementary differentials, which are products of derivatives of the
ODE right-hand side.  Since we haven't specified an ODE, these are indicated
simply by the associated rooted tree.  The rooted trees are printed as nested
lists, essentially in the form used in Butcher's book.  The rooted trees written
in this way can be rendered in LaTex using the package `forest`; unfortunately,
there is no easy way to render them in the browser.

Here is a B-series for a 4th-order method, expanded up to 5th-order terms:


```julia
A = Rational{Int128}[0 0 0 0 0 0 0 0;(-1//6) (1//2) 0 0 0 0 0 0;(-1//10) (1//10) (1//2) 0 0 0 0 0;(-21463//39375) (21017//26250) (-5//9) (1//2) 0 0 0 0;(-59588//54675) (118717//36450) (-4375//2187) 0 (1//2) 0 0 0;(-19993033//9443328) (28508695//3147776) (-13577105//2360832) (-4090625//3147776) (1136025//3147776) (1//2) 0 0;(367020141781//199294617600) (814214904871//22143846400) (-29834937659//1992946176) (-1983358776875//87689631744) (-6702625935//885753856) (688576//109395) (1//2) 0;(1081252805//134140608) (2639189439//74522560) (33646441//4191894) (-7873511875//210792384) (-504040617//14904512) (2110843561//115277085) (13//7) (1//2)];
b = Rational{Int128}[(1081252805//134140608),(2639189439//74522560),(33646441//4191894),(-7873511875//210792384),(-504040617//14904512),(2110843561//115277085),(13//7),(1//2)];
c = Rational{Int128}[0,(1//3),(1//2),(1//5),(2//3),(3//4),(1//4),1];

coeffs4 = bseries(A,b,c,5)
latexify(coeffs4, cdot=false)
```




$F_{f}\mathopen{}\left( \varnothing \right)\mathclose{} + h F_{f}\mathopen{}\left( \rootedtree[] \right)\mathclose{} + \frac{1}{2} h^{2} F_{f}\mathopen{}\left( \rootedtree[[]] \right)\mathclose{} + \frac{1}{6} h^{3} F_{f}\mathopen{}\left( \rootedtree[[[]]] \right)\mathclose{} + \frac{1}{6} h^{3} F_{f}\mathopen{}\left( \rootedtree[[][]] \right)\mathclose{} + \frac{1}{24} h^{4} F_{f}\mathopen{}\left( \rootedtree[[[[]]]] \right)\mathclose{} + \frac{1}{24} h^{4} F_{f}\mathopen{}\left( \rootedtree[[[][]]] \right)\mathclose{} + \frac{1}{8} h^{4} F_{f}\mathopen{}\left( \rootedtree[[[]][]] \right)\mathclose{} + \frac{1}{24} h^{4} F_{f}\mathopen{}\left( \rootedtree[[][][]] \right)\mathclose{} + \frac{1}{48} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[[]]]]] \right)\mathclose{} + \frac{1}{48} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[][]]]] \right)\mathclose{} + \frac{1}{16} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[]][]]] \right)\mathclose{} + \frac{-21246894637}{49670350848} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[]]][]] \right)\mathclose{} + \frac{1}{48} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[][][]]] \right)\mathclose{} + \frac{-722476128287}{1390769823744} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[][]][]] \right)\mathclose{} + \frac{1970748171909370823}{42730909364772720000} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[]][[]]] \right)\mathclose{} + \frac{3898363669}{40242182400} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[]][][]] \right)\mathclose{} + \frac{1}{48} h^{5} F_{f}\mathopen{}\left( \rootedtree[[][][][]] \right)\mathclose{}$



We can also print out the B-series coefficients this way:


```julia
coeffs4
```




    TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Rational{Int128}} with 18 entries:
      RootedTree{Int64}: Int64[]         => 1//1
      RootedTree{Int64}: [1]             => 1//1
      RootedTree{Int64}: [1, 2]          => 1//2
      RootedTree{Int64}: [1, 2, 3]       => 1//6
      RootedTree{Int64}: [1, 2, 2]       => 1//3
      RootedTree{Int64}: [1, 2, 3, 4]    => 1//24
      RootedTree{Int64}: [1, 2, 3, 3]    => 1//12
      RootedTree{Int64}: [1, 2, 3, 2]    => 1//8
      RootedTree{Int64}: [1, 2, 2, 2]    => 1//4
      RootedTree{Int64}: [1, 2, 3, 4, 5] => 1//48
      RootedTree{Int64}: [1, 2, 3, 4, 4] => 1//24
      RootedTree{Int64}: [1, 2, 3, 4, 3] => 1//16
      RootedTree{Int64}: [1, 2, 3, 4, 2] => -21246894637//49670350848
      RootedTree{Int64}: [1, 2, 3, 3, 3] => 1//8
      RootedTree{Int64}: [1, 2, 3, 3, 2] => -722476128287//695384911872
      RootedTree{Int64}: [1, 2, 3, 2, 3] => 1970748171909370823//213654546823863600…
      RootedTree{Int64}: [1, 2, 3, 2, 2] => 3898363669//20121091200
      RootedTree{Int64}: [1, 2, 2, 2, 2] => 1//2



In this form, the rooted trees are printed as level sequences.  The corresponding coefficients are on the right.

# Exact series and local error

We can also get the B-series of the exact solution:


```julia
coeffs_ex = ExactSolution(coeffs4)
```




    TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Rational{Int128}} with 18 entries:
      RootedTree{Int64}: Int64[]         => 1//1
      RootedTree{Int64}: [1]             => 1//1
      RootedTree{Int64}: [1, 2]          => 1//2
      RootedTree{Int64}: [1, 2, 3]       => 1//6
      RootedTree{Int64}: [1, 2, 2]       => 1//3
      RootedTree{Int64}: [1, 2, 3, 4]    => 1//24
      RootedTree{Int64}: [1, 2, 3, 3]    => 1//12
      RootedTree{Int64}: [1, 2, 3, 2]    => 1//8
      RootedTree{Int64}: [1, 2, 2, 2]    => 1//4
      RootedTree{Int64}: [1, 2, 3, 4, 5] => 1//120
      RootedTree{Int64}: [1, 2, 3, 4, 4] => 1//60
      RootedTree{Int64}: [1, 2, 3, 4, 3] => 1//40
      RootedTree{Int64}: [1, 2, 3, 4, 2] => 1//30
      RootedTree{Int64}: [1, 2, 3, 3, 3] => 1//20
      RootedTree{Int64}: [1, 2, 3, 3, 2] => 1//15
      RootedTree{Int64}: [1, 2, 3, 2, 3] => 1//20
      RootedTree{Int64}: [1, 2, 3, 2, 2] => 1//10
      RootedTree{Int64}: [1, 2, 2, 2, 2] => 1//5




```julia
latexify(coeffs_ex,cdot=false)
```




$F_{f}\mathopen{}\left( \varnothing \right)\mathclose{} + h F_{f}\mathopen{}\left( \rootedtree[] \right)\mathclose{} + \frac{1}{2} h^{2} F_{f}\mathopen{}\left( \rootedtree[[]] \right)\mathclose{} + \frac{1}{6} h^{3} F_{f}\mathopen{}\left( \rootedtree[[[]]] \right)\mathclose{} + \frac{1}{6} h^{3} F_{f}\mathopen{}\left( \rootedtree[[][]] \right)\mathclose{} + \frac{1}{24} h^{4} F_{f}\mathopen{}\left( \rootedtree[[[[]]]] \right)\mathclose{} + \frac{1}{24} h^{4} F_{f}\mathopen{}\left( \rootedtree[[[][]]] \right)\mathclose{} + \frac{1}{8} h^{4} F_{f}\mathopen{}\left( \rootedtree[[[]][]] \right)\mathclose{} + \frac{1}{24} h^{4} F_{f}\mathopen{}\left( \rootedtree[[][][]] \right)\mathclose{} + \frac{1}{120} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[[]]]]] \right)\mathclose{} + \frac{1}{120} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[][]]]] \right)\mathclose{} + \frac{1}{40} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[]][]]] \right)\mathclose{} + \frac{1}{30} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[]]][]] \right)\mathclose{} + \frac{1}{120} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[][][]]] \right)\mathclose{} + \frac{1}{30} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[][]][]] \right)\mathclose{} + \frac{1}{40} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[]][[]]] \right)\mathclose{} + \frac{1}{20} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[]][][]] \right)\mathclose{} + \frac{1}{120} h^{5} F_{f}\mathopen{}\left( \rootedtree[[][][][]] \right)\mathclose{}$



We can find the local error by subtracting the exact solution B-series from the RK method B-series:


```julia
latexify(coeffs4-coeffs_ex,cdot=false)
```




$\frac{1}{80} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[[]]]]] \right)\mathclose{} + \frac{1}{80} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[][]]]] \right)\mathclose{} + \frac{3}{80} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[]][]]] \right)\mathclose{} + \frac{-114512864993}{248351754240} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[[]]][]] \right)\mathclose{} + \frac{1}{80} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[][][]]] \right)\mathclose{} + \frac{-3844175612059}{6953849118720} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[][]][]] \right)\mathclose{} + \frac{902475437790052823}{42730909364772720000} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[]][[]]] \right)\mathclose{} + \frac{1886254549}{40242182400} h^{5} F_{f}\mathopen{}\left( \rootedtree[[[]][][]] \right)\mathclose{} + \frac{1}{80} h^{5} F_{f}\mathopen{}\left( \rootedtree[[][][][]] \right)\mathclose{}$



This confirms that the method is of 4th order, since all terms involving
smaller powers of $h$ vanish exactly.  We don't see the $h^6$ and higher
order terms since we only generated the truncated B-series up to 5th order.

For the 2nd-order method, we get:


```julia
latexify(coeffs2-coeffs_ex,cdot=false)
```




$\frac{ - h^{3}}{6} F_{f}\mathopen{}\left( \rootedtree[[[]]] \right)\mathclose{} + h^{3} \left( \frac{-1}{6} + \frac{1}{8 \alpha} \right) F_{f}\mathopen{}\left( \rootedtree[[][]] \right)\mathclose{}$



This confirms again the accuracy of the method, and shows us that we
can eliminate one of the leading error terms completely if we take
$\alpha=3/4$ (this is known as Ralston's method, or sometimes as Heun's method).

# B-Series for a specific ODE

Next, let us define an ODE.  We'll consider the Prothero-Robinson problem:

$$
    y'(t) = \lambda(y-\sin(t)) + \cos(t).
$$

For a non-autonomous ODE like this, it's convenient to rewrite the problem
in autonomous form.  We set $u=[y,t]^T$ and

\begin{align}
u_1'(t) & = \lambda(u_1 - \sin(u_2)) + \cos(u_2) \\
u_2'(t) & = t.
\end{align}


```julia
λ = sp.symbols("λ", real=true)
y, t = sp.symbols("y t", real=true)
h = sp.symbols("h", real=true)

u = [y, t]
ff = [λ*(u[1]-sin(t))+cos(t), 1]
```




$\left[ \begin{array}{r}λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\\1\end{array} \right]$




Finally, we get the B-Series for our RK method applied to our ODE:


```julia
evaluate(ff,u,h,coeffs4)[1]
```




$\frac{h^{5} λ^{3} \left(λ \left(λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\right) - λ \cos{\left(t \right)} - \sin{\left(t \right)}\right)}{48} + \frac{h^{5} λ^{2} \left(λ \sin{\left(t \right)} - \cos{\left(t \right)}\right)}{48} + \frac{h^{5} λ \left(λ \cos{\left(t \right)} + \sin{\left(t \right)}\right)}{48} + \frac{h^{5} \left(- λ \sin{\left(t \right)} + \cos{\left(t \right)}\right)}{48} + \frac{h^{4} λ^{2} \left(λ \left(λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\right) - λ \cos{\left(t \right)} - \sin{\left(t \right)}\right)}{24} + \frac{h^{4} λ \left(λ \sin{\left(t \right)} - \cos{\left(t \right)}\right)}{24} + \frac{h^{4} \left(λ \cos{\left(t \right)} + \sin{\left(t \right)}\right)}{24} + \frac{h^{3} λ \left(λ \left(λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\right) - λ \cos{\left(t \right)} - \sin{\left(t \right)}\right)}{6} + \frac{h^{3} \left(λ \sin{\left(t \right)} - \cos{\left(t \right)}\right)}{6} + \frac{h^{2} \left(λ \left(λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\right) - λ \cos{\left(t \right)} - \sin{\left(t \right)}\right)}{2} + h \left(λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\right) + 1$



Notice that the series is truncated at the same order that we specified
when we initially generated it from the RK coefficients.

Here's the B-Series for the exact solution of the same ODE:


```julia
evaluate(ff,u,h,coeffs_ex)[1]
```




$\frac{h^{5} λ^{3} \left(λ \left(λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\right) - λ \cos{\left(t \right)} - \sin{\left(t \right)}\right)}{120} + \frac{h^{5} λ^{2} \left(λ \sin{\left(t \right)} - \cos{\left(t \right)}\right)}{120} + \frac{h^{5} λ \left(λ \cos{\left(t \right)} + \sin{\left(t \right)}\right)}{120} + \frac{h^{5} \left(- λ \sin{\left(t \right)} + \cos{\left(t \right)}\right)}{120} + \frac{h^{4} λ^{2} \left(λ \left(λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\right) - λ \cos{\left(t \right)} - \sin{\left(t \right)}\right)}{24} + \frac{h^{4} λ \left(λ \sin{\left(t \right)} - \cos{\left(t \right)}\right)}{24} + \frac{h^{4} \left(λ \cos{\left(t \right)} + \sin{\left(t \right)}\right)}{24} + \frac{h^{3} λ \left(λ \left(λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\right) - λ \cos{\left(t \right)} - \sin{\left(t \right)}\right)}{6} + \frac{h^{3} \left(λ \sin{\left(t \right)} - \cos{\left(t \right)}\right)}{6} + \frac{h^{2} \left(λ \left(λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\right) - λ \cos{\left(t \right)} - \sin{\left(t \right)}\right)}{2} + h \left(λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\right) + 1$



And their difference, which is the local error:


```julia
expr = sp.simplify(evaluate(ff,u,h,coeffs4)-evaluate(ff,u,h,coeffs_ex))[1]
```




$\frac{h^{5} λ^{3} \left(λ \left(λ \left(y - \sin{\left(t \right)}\right) + \cos{\left(t \right)}\right) - λ \cos{\left(t \right)} - \sin{\left(t \right)}\right)}{80} + \frac{h^{5} λ^{2} \left(λ \sin{\left(t \right)} - \cos{\left(t \right)}\right)}{80} + \frac{h^{5} λ \left(λ \cos{\left(t \right)} + \sin{\left(t \right)}\right)}{80} + \frac{h^{5} \left(- λ \sin{\left(t \right)} + \cos{\left(t \right)}\right)}{80}$



Because of the simplicity of the PR problem, the error has a relatively simple form:


```julia
sp.simplify(sp.expand(expr))
```




$\frac{h^{5} \left(y λ^{5} - λ^{5} \sin{\left(t \right)} + \cos{\left(t \right)}\right)}{80}$




```julia
sp.collect(sp.expand(expr),λ)
```




$\frac{h^{5} \cos{\left(t \right)}}{80} + λ^{5} \left(\frac{h^{5} y}{80} - \frac{h^{5} \sin{\left(t \right)}}{80}\right)$



# B-series for a generic RK method

We can also examine just the elementary differentials, without specifying a RK method:


```julia
elementary_differentials(ff,u,5)
```




    OrderedDict{RootedTree{Int64, Vector{Int64}}, Vector{SymPy.Sym}} with 18 entries:
      RootedTree{Int64}: Int64[]         => [1, 1]
      RootedTree{Int64}: [1]             => [λ*(y - sin(t)) + cos(t), 1]
      RootedTree{Int64}: [1, 2]          => [λ*(λ*(y - sin(t)) + cos(t)) - λ*cos(t)…
      RootedTree{Int64}: [1, 2, 3]       => [λ*(λ*(λ*(y - sin(t)) + cos(t)) - λ*cos…
      RootedTree{Int64}: [1, 2, 2]       => [λ*sin(t) - cos(t), 0]
      RootedTree{Int64}: [1, 2, 3, 4]    => [λ^2*(λ*(λ*(y - sin(t)) + cos(t)) - λ*c…
      RootedTree{Int64}: [1, 2, 3, 3]    => [λ*(λ*sin(t) - cos(t)), 0]
      RootedTree{Int64}: [1, 2, 3, 2]    => [0, 0]
      RootedTree{Int64}: [1, 2, 2, 2]    => [λ*cos(t) + sin(t), 0]
      RootedTree{Int64}: [1, 2, 3, 4, 5] => [λ^3*(λ*(λ*(y - sin(t)) + cos(t)) - λ*c…
      RootedTree{Int64}: [1, 2, 3, 4, 4] => [λ^2*(λ*sin(t) - cos(t)), 0]
      RootedTree{Int64}: [1, 2, 3, 4, 3] => [0, 0]
      RootedTree{Int64}: [1, 2, 3, 4, 2] => [0, 0]
      RootedTree{Int64}: [1, 2, 3, 3, 3] => [λ*(λ*cos(t) + sin(t)), 0]
      RootedTree{Int64}: [1, 2, 3, 3, 2] => [0, 0]
      RootedTree{Int64}: [1, 2, 3, 2, 3] => [0, 0]
      RootedTree{Int64}: [1, 2, 3, 2, 2] => [0, 0]
      RootedTree{Int64}: [1, 2, 2, 2, 2] => [-λ*sin(t) + cos(t), 0]


