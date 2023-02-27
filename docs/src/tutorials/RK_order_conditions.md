# Runge-Kutta order conditions

In this tutorial, we generate order conditions for a generic explicit Runge-Kutta method.

First, we create symbolic coefficient arrays with the appropriate structure.
There are several symbolic packages you can use. Here, we will use SymPy.jl.


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


The output above is truncated.  We can see the full expressions as follows:


```@example bseries-RK-order-conditions
for (tree, order_condition) in error
    println(order_condition)
end
```

