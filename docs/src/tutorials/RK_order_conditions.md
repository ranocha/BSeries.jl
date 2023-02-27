# Runge-Kutta order conditions

In this tutorial, we generate order conditions for a generic explicit Runge-Kutta method.

First, we create symbolic coefficient arrays with the appropriate structure.


```@example bseries-RK-order-conditions
using BSeries, SymPy, Latexify

s = 4  # Stages
p = 4  # Desired order of accuracy

A = Array{Sym,2}(undef, s, s)
b = Array{Sym,1}(undef, s)
c = Array{Sym,1}(undef, s);
```


```@example bseries-RK-order-conditions
for i in 1:s
    b[i] = symbols("b$i", real = true)
    for j in 1:i-1
        A[i, j] = symbols("a$i$j", real = true)
    end
    for j in i:s
        A[i, j] = 0
    end
end

for i in 1:s
    c[i] = 0
    for j in 1:i-1
        c[i] += A[i, j]
    end
end
```

Next we generate the B-series for the RK method and the exact solution, and take their difference.


```@example bseries-RK-order-conditions

rk3 = bseries(A, b, c, p)
exact = ExactSolution(rk3)
error = rk3 - exact
```




    TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Sym} with 9 entries:
      ∅       => 0
      τ       => b1 + b2 + b3 + b4 - 1
      [τ]     => a21*b2 + b3*(a31 + a32) + b4*(a41 + a42 + a43) - 1/2
      [[τ]]   => a21*a32*b3 + b4*(a21*a42 + a43*(a31 + a32)) - 1/6
      [τ²]    => a21^2*b2 + b3*(a31 + a32)^2 + b4*(a41 + a42 + a43)^2 - 1/3
      [[[τ]]] => a21*a32*a43*b4 - 1/24
      [[τ²]]  => a21^2*a32*b3 + b4*(a21^2*a42 + a43*(a31 + a32)^2) - 1/12
      [[τ]τ]  => a21*a32*b3*(a31 + a32) + b4*(a21*a42 + a43*(a31 + a32))*(a41 + a42…
      [τ³]    => a21^3*b2 + b3*(a31 + a32)^3 + b4*(a41 + a42 + a43)^3 - 1/4



The output above is truncated in the Jupyter notebook.  We can see the full expressions as follows:


```julia
for (tree,order_condition) in error
    println(order_condition)
end
```

    0
    b1 + b2 + b3 + b4 - 1
    a21*b2 + b3*(a31 + a32) + b4*(a41 + a42 + a43) - 1/2
    a21*a32*b3 + b4*(a21*a42 + a43*(a31 + a32)) - 1/6
    a21^2*b2 + b3*(a31 + a32)^2 + b4*(a41 + a42 + a43)^2 - 1/3
    a21*a32*a43*b4 - 1/24
    a21^2*a32*b3 + b4*(a21^2*a42 + a43*(a31 + a32)^2) - 1/12
    a21*a32*b3*(a31 + a32) + b4*(a21*a42 + a43*(a31 + a32))*(a41 + a42 + a43) - 1/8
    a21^3*b2 + b3*(a31 + a32)^3 + b4*(a41 + a42 + a43)^3 - 1/4

