using Test
using BSeries

using Latexify: latexify

using StaticArrays: @SArray

using SymEngine: SymEngine
using SymPy: SymPy
using Symbolics: Symbolics


@testset "BSeries" begin


@testset "lazy representation of exact ODE solution" begin
  exact = ExactSolution{Rational{Int}}()
  terms = collect(Iterators.take(exact, 4))
  @test terms == [1//1, 1//2, 1//6, 1//3]
  @test exact == ExactSolution(exact)
end


@testset "latexify" begin
  # explicit midpoint method
  A = @SArray [0 0; 1//2 0]; b = @SArray [0, 1//1]; c = @SArray [0, 1//2];

  series_integrator = bseries(A, b, c, 2)
  @test_nowarn latexify(series_integrator)
  @test_nowarn latexify(series_integrator, cdot=false)
  @test_nowarn latexify(series_integrator, dt=SymEngine.symbols("h"))
end


@testset "AbstractDict interface" begin
  # explicit midpoint method
  A = @SArray [0 0; 1//2 0]; b = @SArray [0, 1//1]; c = @SArray [0, 1//2];
  series_integrator = bseries(A, b, c, 2)

  # These are mostly simple smoke tests
  @test_nowarn begin
    @test isempty(empty(series_integrator, RootedTrees.RootedTree{Int32, Vector{Int32}}, Rational{Int32}))
  end
  empty!(series_integrator)
  @test isempty(series_integrator)
  @test isempty(keys(series_integrator))
  @test isempty(values(series_integrator))
  @test length(series_integrator) == 0
  @test get(series_integrator, rootedtree([1]), :default) == :default
  @test get(series_integrator, rootedtree([1])) do
    :default
  end == :default
  @test getkey(series_integrator, rootedtree([1]), :default) == :default
  @test_nowarn sizehint!(series_integrator, 10)

  @test get!(series_integrator, rootedtree([1]), 2) == 2
  @test pop!(series_integrator, rootedtree([1])) == 2
  @test pop!(series_integrator, rootedtree([1]), :default) == :default
  series_integrator[rootedtree([1])] = 3
  @test isempty(delete!(series_integrator, rootedtree([1])))
end


@testset "substitute" begin
  Symbolics.@variables a1 a2 a31 a32
  a = OrderedDict{RootedTrees.RootedTree{Int, Vector{Int}}, Symbolics.Num}(
        rootedtree(Int[])     => 1,
        rootedtree([1])       => a1,
        rootedtree([1, 2])    => a2,
        rootedtree([1, 2, 2]) => a31,
        rootedtree([1, 2, 3]) => a32)

  Symbolics.@variables b1 b2 b31 b32
  b = OrderedDict{RootedTrees.RootedTree{Int, Vector{Int}}, Symbolics.Num}(
        rootedtree(Int[])     => 0,
        rootedtree([1])       => b1,
        rootedtree([1, 2])    => b2,
        rootedtree([1, 2, 2]) => b31,
        rootedtree([1, 2, 3]) => b32)

  # See equation (6) of
  # - Philippe Chartier, Ernst Hairer and Gilles Vilmart (2007)
  #   Numerical integrators based on modified differential equations
  #   [DOI: 10.1090/S0025-5718-07-01967-9](https://doi.org/10.1090/S0025-5718-07-01967-9)

  t = rootedtree([1])
  @inferred substitute(b, a, t)
  @test isequal(substitute(b, a, t), a1 * b1)

  t = rootedtree([1, 2])
  @inferred substitute(b, a, t)
  @test isequal(substitute(b, a, t), a1 * b2 + a2 * b1^2)

  t = rootedtree([1, 2, 2])
  @inferred substitute(b, a, t)
  @test isequal(substitute(b, a, t), a1 * b31 + 2 * a2 * b1 * b2 + a31 * b1^3)

  t = rootedtree([1, 2, 3])
  @inferred substitute(b, a, t)
  @test isequal(substitute(b, a, t), a1 * b32 + 2 * a2 * b1 * b2 + a32 * b1^3)
end


@testset "compose" begin
  Symbolics.@variables a0 a1 a2 a31 a32
  a = OrderedDict{RootedTrees.RootedTree{Int, Vector{Int}}, Symbolics.Num}(
        rootedtree(Int[])     => a0,
        rootedtree([1])       => a1,
        rootedtree([1, 2])    => a2,
        rootedtree([1, 2, 2]) => a31,
        rootedtree([1, 2, 3]) => a32)

  Symbolics.@variables b1 b2 b31 b32
  b = OrderedDict{RootedTrees.RootedTree{Int, Vector{Int}}, Symbolics.Num}(
        rootedtree(Int[])     => 0,
        rootedtree([1])       => b1,
        rootedtree([1, 2])    => b2,
        rootedtree([1, 2, 2]) => b31,
        rootedtree([1, 2, 3]) => b32)

  # See equation (5) of
  # - Philippe Chartier, Ernst Hairer and Gilles Vilmart (2007)
  #   Numerical integrators based on modified differential equations
  #   [DOI: 10.1090/S0025-5718-07-01967-9](https://doi.org/10.1090/S0025-5718-07-01967-9)
  t = rootedtree([1])
  @inferred compose(b, a, t)
  @test isequal(compose(b, a, t), a0 * b1 + a1)

  t = rootedtree([1, 2])
  @inferred compose(b, a, t)
  @test isequal(compose(b, a, t), a0 * b2 + a1 * b1 + a2)

  t = rootedtree([1, 2, 2])
  @inferred compose(b, a, t)
  @test isequal(compose(b, a, t), a0 * b31 + a1 * b1^2 + 2 * a2 * b1 + a31)

  t = rootedtree([1, 2, 3])
  @inferred compose(b, a, t)
  @test isequal(compose(b, a, t), a0 * b32 + a1 * b2 + a2 * b1 + a32)
end


@testset "modified_equation" begin
  t1  = rootedtree([1])
  t2  = rootedtree([1, 2])
  t31 = rootedtree([1, 2, 2])
  t32 = rootedtree([1, 2, 3])
  t41 = rootedtree([1, 2, 2, 2])
  t42 = rootedtree([1, 2, 2, 3])
  t43 = rootedtree([1, 2, 3, 3])
  t44 = rootedtree([1, 2, 3, 4])
  t51 = rootedtree([1, 2, 2, 2, 2])
  t52 = rootedtree([1, 2, 2, 2, 3])
  t53 = rootedtree([1, 2, 2, 3, 3])
  t54 = rootedtree([1, 2, 2, 3, 4])
  t55 = rootedtree([1, 2, 3, 2, 3])
  t56 = rootedtree([1, 2, 3, 3, 3])
  t57 = rootedtree([1, 2, 3, 3, 4])
  t58 = rootedtree([1, 2, 3, 4, 4])
  t59 = rootedtree([1, 2, 3, 4, 5])

  @testset "SSPRK(2, 2) aka Heun's method" begin
    A = @SArray [0.0 0.0;
                1.0 0.0]
    b = @SArray [1/2, 1/2]
    c = @SArray [0.0, 1.0]
    order = 5
    series = modified_equation(A, b, c, order)

    # tested with the Python package BSeries
    @test series[t1 ] ≈ 1.0                   atol=10*eps()
    @test series[t2 ] ≈ 0                     atol=10*eps()
    @test series[t31] ≈ 0.166666666666667     atol=10*eps()
    @test series[t32] ≈ -0.166666666666667    atol=10*eps()
    @test series[t41] ≈ -1.11022302462516e-16 atol=10*eps()
    @test series[t42] ≈ -0.125000000000000    atol=10*eps()
    @test series[t43] ≈ -2.77555756156289e-17 atol=10*eps()
    @test series[t44] ≈ 0.125000000000000     atol=10*eps()
    @test series[t51] ≈ -0.0333333333333332   atol=10*eps()
    @test series[t52] ≈ -0.0583333333333333   atol=10*eps()
    @test series[t55] ≈ 0.0750000000000000    atol=10*eps()
    @test series[t53] ≈ 0.0583333333333333    atol=10*eps()
    @test series[t56] ≈ 0.0333333333333334    atol=10*eps()
    @test series[t54] ≈ 0.0500000000000000    atol=10*eps()
    @test series[t57] ≈ 0.0583333333333333    atol=10*eps()
    @test series[t58] ≈ -0.0583333333333333   atol=10*eps()
    @test series[t59] ≈ -0.0500000000000000   atol=10*eps()
  end


  @testset "SSPRK(3, 3)" begin
    A = @SArray [0.0 0.0 0.0;
                1.0 0.0 0.0;
                1/4 1/4 0.0]
    b = @SArray [1/6, 1/6, 2/3]
    c = @SArray [0.0, 1.0, 1/2]
    order = 5
    series = modified_equation(A, b, c, order)

    # tested with the Python package BSeries
    @test series[t1 ] ≈ 1.0                  atol=10*eps()
    @test series[t2 ] ≈ 0                    atol=10*eps()
    @test series[t31] ≈ 0                    atol=10*eps()
    @test series[t32] ≈ 0                    atol=10*eps()
    @test series[t41] ≈ 0                    atol=10*eps()
    @test series[t42] ≈ -0.0416666666666667  atol=10*eps()
    @test series[t43] ≈ 0.0833333333333333   atol=10*eps()
    @test series[t44] ≈ -0.0416666666666667  atol=10*eps()
    @test series[t51] ≈ 0.00833333333333330  atol=10*eps()
    @test series[t52] ≈ -0.0166666666666667  atol=10*eps()
    @test series[t55] ≈ 0.0333333333333333   atol=10*eps()
    @test series[t53] ≈ 0.0166666666666667   atol=10*eps()
    @test series[t56] ≈ -0.00833333333333333 atol=10*eps()
    @test series[t54] ≈ 0.00833333333333333  atol=10*eps()
    @test series[t57] ≈ -0.0250000000000000  atol=10*eps()
    @test series[t58] ≈ -0.0166666666666667  atol=10*eps()
    @test series[t59] ≈ 0.0333333333333333   atol=10*eps()
  end


  @testset "Explicit midpoint method (Runge's method)" begin
    A = @SArray [0 0; 1//2 0]
    b = @SArray [0, 1//1]
    c = @SArray [0, 1//2]
    order = 5
    series = modified_equation(A, b, c, order)

    # tested with the Python package BSeries
    @test series[t1 ] == 1
    @test series[t2 ] == 0
    @test series[t31] == -1//12
    @test series[t32] == -1//6
    @test series[t41] == 0
    @test series[t42] == 0
    @test series[t43] == 1//8
    @test series[t44] == 1//8
    @test series[t51] == 7//240
    @test series[t52] == 1//40
    @test series[t55] == 1//30
    @test series[t53] == 3//80
    @test series[t56] == -7//240
    @test series[t54] == 7//240
    @test series[t57] == -1//40
    @test series[t58] == -19//240
    @test series[t59] == -1//20
  end
end


@testset "modified_equation with elementary differentials" begin
  @testset "SymEngine.jl" begin
    # Lotka-Volterra model
    dt = SymEngine.symbols("dt")
    p, q = u = SymEngine.symbols("p, q")
    f = [p * (2 - q), q * (p - 1)]

    # Explicit Euler method
    A = @SArray [0//1;]
    b = @SArray [1//1]
    c = @SArray [0//1]

    # tested with the Python package BSeries
    series2 = modified_equation(f, u, dt, A, b, c, 2)
    series2_reference = [
      -dt*(-p*q*(p - 1) + p*(2 - q)^2)/2 + p*(2 - q),
      -dt*( p*q*(2 - q) + q*(p - 1)^2)/2 + q*(p - 1)
    ]
    @test mapreduce(isequal, &, series2, series2_reference)

    # tested with the Python package BSeries
    series3 = modified_equation(f, u, dt, A, b, c, 3)
    series3_reference = [
      -dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p^2*q*(2 - q) - p*q*(2 - q)*(p - 1) - p*q*(p - 1)^2 + p*(2 - q)^3)/3 - dt*(-p*q*(p - 1) + p*(2 - q)^2)/2 + p*(2 - q),
      dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p*q^2*(p - 1) + p*q*(2 - q)^2 + p*q*(2 - q)*(p - 1) + q*(p - 1)^3)/3 - dt*(p*q*(2 - q) + q*(p - 1)^2)/2 + q*(p - 1)
    ]
    @test mapreduce(iszero ∘ SymEngine.expand, &, series3 - series3_reference)

    # tested with the Python package BSeries
    series4 = modified_equation(f, u, dt, A, b, c, 4)
    series4_reference = [
      -dt^3*(-2*p^2*q*(2 - q)*(p - 1) + 2*p*q*(2 - q)*(p - 1)*(q - 2))/12 - dt^3*(-p^2*q*(2 - q)^2 + p*q^2*(p - 1)^2 + p*q*(1 - p)*(2 - q)*(p - 1) + p*q*(2 - q)*(p - 1)*(q - 2))/12 - dt^3*(p^2*q^2*(p - 1) - 2*p^2*q*(2 - q)^2 - p^2*q*(2 - q)*(p - 1) - p*q*(2 - q)^2*(p - 1) - p*q*(2 - q)*(p - 1)^2 - p*q*(p - 1)^3 + p*(2 - q)^4)/4 - dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p^2*q*(2 - q) - p*q*(2 - q)*(p - 1) - p*q*(p - 1)^2 + p*(2 - q)^3)/3 - dt*(-p*q*(p - 1) + p*(2 - q)^2)/2 + p*(2 - q),
        -dt^3*(-2*p*q^2*(2 - q)*(p - 1) + 2*p*q*(2 - q)*(p - 1)^2)/12 - dt^3*(p^2*q*(2 - q)^2 - p*q^2*(p - 1)^2 + p*q*(2 - q)^2*(p - 1) + p*q*(2 - q)*(p - 1)^2)/12 - dt^3*(-p^2*q^2*(2 - q) - p*q^2*(2 - q)*(p - 1) - 2*p*q^2*(p - 1)^2 + p*q*(2 - q)^3 + p*q*(2 - q)^2*(p - 1) + p*q*(2 - q)*(p - 1)^2 + q*(p - 1)^4)/4 + dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p*q^2*(p - 1) + p*q*(2 - q)^2 + p*q*(2 - q)*(p - 1) + q*(p - 1)^3)/3 - dt*(p*q*(2 - q) + q*(p - 1)^2)/2 + q*(p - 1)
    ]
    @test mapreduce(iszero ∘ SymEngine.expand, &, series4 - series4_reference)
  end

  @testset "SymPy.jl" begin
    # Lotka-Volterra model
    dt = SymPy.symbols("dt")
    p, q = u = SymPy.symbols("p, q")
    f = [p * (2 - q), q * (p - 1)]

    # Explicit Euler method
    A = @SArray [0//1;]
    b = @SArray [1//1]
    c = @SArray [0//1]

    # tested with the Python package BSeries
    series2 = modified_equation(f, u, dt, A, b, c, 2)
    series2_reference = [
      -dt*(-p*q*(p - 1) + p*(2 - q)^2)/2 + p*(2 - q),
      -dt*( p*q*(2 - q) + q*(p - 1)^2)/2 + q*(p - 1)
    ]
    @test mapreduce(isequal, &, series2, series2_reference)

    # tested with the Python package BSeries
    series3 = modified_equation(f, u, dt, A, b, c, 3)
    series3_reference = [
      -dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p^2*q*(2 - q) - p*q*(2 - q)*(p - 1) - p*q*(p - 1)^2 + p*(2 - q)^3)/3 - dt*(-p*q*(p - 1) + p*(2 - q)^2)/2 + p*(2 - q),
      dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p*q^2*(p - 1) + p*q*(2 - q)^2 + p*q*(2 - q)*(p - 1) + q*(p - 1)^3)/3 - dt*(p*q*(2 - q) + q*(p - 1)^2)/2 + q*(p - 1)
    ]
    @test mapreduce(iszero ∘ SymPy.expand, &, series3 - series3_reference)

    # tested with the Python package BSeries
    series4 = modified_equation(f, u, dt, A, b, c, 4)
    series4_reference = [
      -dt^3*(-2*p^2*q*(2 - q)*(p - 1) + 2*p*q*(2 - q)*(p - 1)*(q - 2))/12 - dt^3*(-p^2*q*(2 - q)^2 + p*q^2*(p - 1)^2 + p*q*(1 - p)*(2 - q)*(p - 1) + p*q*(2 - q)*(p - 1)*(q - 2))/12 - dt^3*(p^2*q^2*(p - 1) - 2*p^2*q*(2 - q)^2 - p^2*q*(2 - q)*(p - 1) - p*q*(2 - q)^2*(p - 1) - p*q*(2 - q)*(p - 1)^2 - p*q*(p - 1)^3 + p*(2 - q)^4)/4 - dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p^2*q*(2 - q) - p*q*(2 - q)*(p - 1) - p*q*(p - 1)^2 + p*(2 - q)^3)/3 - dt*(-p*q*(p - 1) + p*(2 - q)^2)/2 + p*(2 - q),
        -dt^3*(-2*p*q^2*(2 - q)*(p - 1) + 2*p*q*(2 - q)*(p - 1)^2)/12 - dt^3*(p^2*q*(2 - q)^2 - p*q^2*(p - 1)^2 + p*q*(2 - q)^2*(p - 1) + p*q*(2 - q)*(p - 1)^2)/12 - dt^3*(-p^2*q^2*(2 - q) - p*q^2*(2 - q)*(p - 1) - 2*p*q^2*(p - 1)^2 + p*q*(2 - q)^3 + p*q*(2 - q)^2*(p - 1) + p*q*(2 - q)*(p - 1)^2 + q*(p - 1)^4)/4 + dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p*q^2*(p - 1) + p*q*(2 - q)^2 + p*q*(2 - q)*(p - 1) + q*(p - 1)^3)/3 - dt*(p*q*(2 - q) + q*(p - 1)^2)/2 + q*(p - 1)
    ]
    @test mapreduce(iszero ∘ SymPy.expand, &, series4 - series4_reference)
  end

  @testset "Symbolics.jl" begin
    # Lotka-Volterra model
    Symbolics.@variables dt
    u = Symbolics.@variables p q
    f = [p * (2 - q), q * (p - 1)]

    # Explicit Euler method
    A = @SArray [0//1;]
    b = @SArray [1//1]
    c = @SArray [0//1]

    # tested with the Python package BSeries
    series2 = modified_equation(f, u, dt, A, b, c, 2)
    series2_reference = [
      -dt*(-p*q*(p - 1) + p*(2 - q)^2)/2 + p*(2 - q),
      -dt*( p*q*(2 - q) + q*(p - 1)^2)/2 + q*(p - 1)
    ]
    @test mapreduce(isequal, &, series2, series2_reference)

    # tested with the Python package BSeries
    series3 = modified_equation(f, u, dt, A, b, c, 3)
    series3_reference = [
      -dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p^2*q*(2 - q) - p*q*(2 - q)*(p - 1) - p*q*(p - 1)^2 + p*(2 - q)^3)/3 - dt*(-p*q*(p - 1) + p*(2 - q)^2)/2 + p*(2 - q),
      dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p*q^2*(p - 1) + p*q*(2 - q)^2 + p*q*(2 - q)*(p - 1) + q*(p - 1)^3)/3 - dt*(p*q*(2 - q) + q*(p - 1)^2)/2 + q*(p - 1)
    ]
    @test mapreduce(iszero ∘ Symbolics.expand, &, series3 - series3_reference)

    # tested with the Python package BSeries
    series4 = modified_equation(f, u, dt, A, b, c, 4)
    series4_reference = [
      -dt^3*(-2*p^2*q*(2 - q)*(p - 1) + 2*p*q*(2 - q)*(p - 1)*(q - 2))/12 - dt^3*(-p^2*q*(2 - q)^2 + p*q^2*(p - 1)^2 + p*q*(1 - p)*(2 - q)*(p - 1) + p*q*(2 - q)*(p - 1)*(q - 2))/12 - dt^3*(p^2*q^2*(p - 1) - 2*p^2*q*(2 - q)^2 - p^2*q*(2 - q)*(p - 1) - p*q*(2 - q)^2*(p - 1) - p*q*(2 - q)*(p - 1)^2 - p*q*(p - 1)^3 + p*(2 - q)^4)/4 - dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p^2*q*(2 - q) - p*q*(2 - q)*(p - 1) - p*q*(p - 1)^2 + p*(2 - q)^3)/3 - dt*(-p*q*(p - 1) + p*(2 - q)^2)/2 + p*(2 - q),
        -dt^3*(-2*p*q^2*(2 - q)*(p - 1) + 2*p*q*(2 - q)*(p - 1)^2)/12 - dt^3*(p^2*q*(2 - q)^2 - p*q^2*(p - 1)^2 + p*q*(2 - q)^2*(p - 1) + p*q*(2 - q)*(p - 1)^2)/12 - dt^3*(-p^2*q^2*(2 - q) - p*q^2*(2 - q)*(p - 1) - 2*p*q^2*(p - 1)^2 + p*q*(2 - q)^3 + p*q*(2 - q)^2*(p - 1) + p*q*(2 - q)*(p - 1)^2 + q*(p - 1)^4)/4 + dt^2*p*q*(2 - q)*(p - 1)/6 + dt^2*(-p*q^2*(p - 1) + p*q*(2 - q)^2 + p*q*(2 - q)*(p - 1) + q*(p - 1)^3)/3 - dt*(p*q*(2 - q) + q*(p - 1)^2)/2 + q*(p - 1)
    ]
    @test mapreduce(iszero ∘ Symbolics.expand, &, series4 - series4_reference)
  end
end


@testset "modifying_integrator" begin
  t1  = rootedtree([1])
  t2  = rootedtree([1, 2])
  t31 = rootedtree([1, 2, 2])
  t32 = rootedtree([1, 2, 3])
  t41 = rootedtree([1, 2, 2, 2])
  t42 = rootedtree([1, 2, 2, 3])
  t43 = rootedtree([1, 2, 3, 3])
  t44 = rootedtree([1, 2, 3, 4])
  t51 = rootedtree([1, 2, 2, 2, 2])
  t52 = rootedtree([1, 2, 2, 2, 3])
  t53 = rootedtree([1, 2, 2, 3, 3])
  t54 = rootedtree([1, 2, 2, 3, 4])
  t55 = rootedtree([1, 2, 3, 2, 3])
  t56 = rootedtree([1, 2, 3, 3, 3])
  t57 = rootedtree([1, 2, 3, 3, 4])
  t58 = rootedtree([1, 2, 3, 4, 4])
  t59 = rootedtree([1, 2, 3, 4, 5])

  @testset "SSPRK(2, 2) aka Heun's method" begin
    A = @SArray [0.0 0.0;
                1.0 0.0]
    b = @SArray [1/2, 1/2]
    c = @SArray [0.0, 1.0]
    order = 5
    series = modifying_integrator(A, b, c, order)

    # tested with the Python package BSeries
    @test series[t1 ] ≈ 1                    atol=10*eps()
    @test series[t2 ] ≈ 0                    atol=10*eps()
    @test series[t31] ≈ -0.166666666666667   atol=10*eps()
    @test series[t32] ≈ 1/6                  atol=10*eps()
    @test series[t41] ≈ 8.32667268468867e-17 atol=10*eps()
    @test series[t42] ≈ 0.125000000000000    atol=10*eps()
    @test series[t43] ≈ 1.38777878078145e-17 atol=10*eps()
    @test series[t44] ≈ -0.125000000000000   atol=10*eps()
    @test series[t51] ≈ 0.200000000000000    atol=10*eps()
    @test series[t52] ≈ 0.0583333333333333   atol=10*eps()
    @test series[t55] ≈ 0.00833333333333335  atol=10*eps()
    @test series[t53] ≈ -0.0583333333333333  atol=10*eps()
    @test series[t56] ≈ -0.200000000000000   atol=10*eps()
    @test series[t54] ≈ -0.133333333333333   atol=10*eps()
    @test series[t57] ≈ -0.0583333333333333  atol=10*eps()
    @test series[t58] ≈ 0.0583333333333333   atol=10*eps()
    @test series[t59] ≈ 0.133333333333333    atol=10*eps()
  end


  @testset "SSPRK(3, 3)" begin
    A = @SArray [0.0 0.0 0.0;
                1.0 0.0 0.0;
                1/4 1/4 0.0]
    b = @SArray [1/6, 1/6, 2/3]
    c = @SArray [0.0, 1.0, 1/2]
    order = 5
    series = modifying_integrator(A, b, c, order)

    # tested with the Python package BSeries
    @test series[t1 ] ≈ 1                    atol=10*eps()
    @test series[t2 ] ≈ 0                    atol=10*eps()
    @test series[t31] ≈ 0                    atol=10*eps()
    @test series[t32] ≈ 0                    atol=10*eps()
    @test series[t41] ≈ 0                    atol=10*eps()
    @test series[t42] ≈ 0.0416666666666667   atol=10*eps()
    @test series[t43] ≈ -0.0833333333333333  atol=10*eps()
    @test series[t44] ≈ 1/24                 atol=10*eps()
    @test series[t51] ≈ -0.00833333333333330 atol=10*eps()
    @test series[t52] ≈ 0.0166666666666667   atol=10*eps()
    @test series[t55] ≈ -0.0333333333333333  atol=10*eps()
    @test series[t53] ≈ -0.0166666666666667  atol=10*eps()
    @test series[t56] ≈ 0.00833333333333332  atol=10*eps()
    @test series[t54] ≈ -0.00833333333333334 atol=10*eps()
    @test series[t57] ≈ 0.0250000000000000   atol=10*eps()
    @test series[t58] ≈ 1/60                 atol=10*eps()
    @test series[t59] ≈ -0.0333333333333333  atol=10*eps()
  end


  @testset "Explicit midpoint method (Runge's method)" begin
    A = @SArray [0 0; 1//2 0]
    b = @SArray [0, 1//1]
    c = @SArray [0, 1//2]
    order = 5
    series = modifying_integrator(A, b, c, order)

    # tested with the Python package BSeries
    @test series[t1 ] == 1
    @test series[t2 ] == 0
    @test series[t31] == 1//12
    @test series[t32] == 1//6
    @test series[t41] == 0
    @test series[t42] == 0
    @test series[t43] == -1//8
    @test series[t44] == -1//8
    @test series[t51] == 1//80
    @test series[t52] == 1//60
    @test series[t55] == 7//240
    @test series[t53] == 1//240
    @test series[t56] == 9//80
    @test series[t54] == 1//80
    @test series[t57] == 13//120
    @test series[t58] == 13//80
    @test series[t59] == 2//15
  end
end


@testset "modifying_integrator with elementary differentials" begin
  @testset "SymEngine.jl" begin
    dt = SymEngine.symbols("dt")
    α, β, γ = SymEngine.symbols("alpha, beta, gamma")
    u1, u2, u3 = u = SymEngine.symbols("u1 u2 u3")
    f = [α*u[2]*u[3], β*u[3]*u[1], γ*u[1]*u[2]]

    # Implicit midpoint method
    A = @SArray [1//2;]
    b = @SArray [1//1]
    c = @SArray [1//2]

    # See eq. (12) of
    # - Philippe Chartier, Ernst Hairer and Gilles Vilmart (2007)
    #   Numerical integrators based on modified differential equations
    #   [DOI: 10.1090/S0025-5718-07-01967-9](https://doi.org/10.1090/S0025-5718-07-01967-9)

    series = modifying_integrator(f, u, dt, A, b, c, 5)
    # Dirty workaround since there is no way to get the polynomial coefficients
    # at the moment.
    # We differentiate the expression and set `dt` to zero to get the corresponding
    # coefficient divided by factorial(how often we needed to differentiate).
    s3_reference = -(α * β * u3^2 + α * γ * u2^2 + β * γ * u1^2) / 12
    for i in eachindex(f)
      s3 = SymEngine.subs(
        1//2 * SymEngine.diff(SymEngine.diff((series[i] - f[i]) / f[i], dt), dt),
        dt => 0)
      @test iszero(SymEngine.expand(s3 - s3_reference))
    end

    s5_reference = 6//5 * s3_reference^2 + 1//60 * α * β * γ * (
      β * u1^2 * u3^2 + γ * u2^2 * u1^2 + α * u3^2 * u2^2
    )
    for i in eachindex(f)
      s5 = SymEngine.subs(
        1//24 * SymEngine.diff(SymEngine.diff(SymEngine.diff(SymEngine.diff(
          (series[i] - f[i]) / f[i], dt), dt), dt), dt),
        dt => 0)
      @test iszero(SymEngine.expand(s5 - s5_reference))
    end
  end

  @testset "SymPy.jl" begin
    dt = SymPy.symbols("dt")
    α, β, γ = SymPy.symbols("alpha, beta, gamma")
    u1, u2, u3 = u = SymPy.symbols("u1 u2 u3")
    f = [α*u[2]*u[3], β*u[3]*u[1], γ*u[1]*u[2]]

    # Implicit midpoint method
    A = @SArray [1//2;]
    b = @SArray [1//1]
    c = @SArray [1//2]

    # See eq. (12) of
    # - Philippe Chartier, Ernst Hairer and Gilles Vilmart (2007)
    #   Numerical integrators based on modified differential equations
    #   [DOI: 10.1090/S0025-5718-07-01967-9](https://doi.org/10.1090/S0025-5718-07-01967-9)

    series = modifying_integrator(f, u, dt, A, b, c, 5)
    # Dirty workaround used also for the other symbolic setups - just make
    # it consistent here, although we could use another approach with SymPy.jl.
    # We differentiate the expression and set `dt` to zero to get the corresponding
    # coefficient divided by factorial(how often we needed to differentiate).
    s3_reference = -(α * β * u3^2 + α * γ * u2^2 + β * γ * u1^2) / 12
    for i in eachindex(f)
      s3 = SymPy.subs(
        1//2 * SymPy.diff(SymPy.diff((series[i] - f[i]) / f[i], dt), dt),
        dt => 0)
      @test iszero(SymPy.expand(s3 - s3_reference))
    end

    s5_reference = 6//5 * s3_reference^2 + 1//60 * α * β * γ * (
      β * u1^2 * u3^2 + γ * u2^2 * u1^2 + α * u3^2 * u2^2
    )
    for i in eachindex(f)
      s5 = SymPy.subs(
        1//24 * SymPy.diff(SymPy.diff(SymPy.diff(SymPy.diff(
          (series[i] - f[i]) / f[i], dt), dt), dt), dt),
        dt => 0)
      @test iszero(SymPy.expand(s5 - s5_reference))
    end
  end

  @testset "Symbolics.jl" begin
    Symbolics.@variables dt
    Symbolics.@variables α β γ
    u = Symbolics.@variables u1 u2 u3
    f = [α*u[2]*u[3], β*u[3]*u[1], γ*u[1]*u[2]]

    # Implicit midpoint method
    A = @SArray [1//2;]
    b = @SArray [1//1]
    c = @SArray [1//2]

    # See eq. (12) of
    # - Philippe Chartier, Ernst Hairer and Gilles Vilmart (2007)
    #   Numerical integrators based on modified differential equations
    #   [DOI: 10.1090/S0025-5718-07-01967-9](https://doi.org/10.1090/S0025-5718-07-01967-9)

    series = modifying_integrator(f, u, dt, A, b, c, 5)
    # Dirty workaround since there is no way to get the polynomial coefficients
    # in Symbolics.jl at the moment, see
    # https://github.com/JuliaSymbolics/Symbolics.jl/issues/216
    # We differentiate the expression and set `dt` to zero to get the corresponding
    # coefficient divided by factorial(how often we needed to differentiate).
    d = Symbolics.Differential(dt)
    s3_reference = -(α * β * u3^2 + α * γ * u2^2 + β * γ * u1^2) / 12
    for i in eachindex(f)
      s3 = Symbolics.simplify_fractions(Symbolics.substitute(
        Symbolics.expand_derivatives(1//2 * d(d((series[i] - f[i]) / f[i]))),
        dt => 0))
      @test isequal(s3, s3_reference)
    end

    s5_reference = 6//5 * s3_reference^2 + 1//60 * α * β * γ * (
      β * u1^2 * u3^2 + γ * u2^2 * u1^2 + α * u3^2 * u2^2
    )
    for i in eachindex(f)
      s5 = Symbolics.simplify_fractions(Symbolics.substitute(
        Symbolics.expand_derivatives(1//24 * d(d(d(d((series[i] - f[i]) / f[i]))))),
        dt => 0))
      @test iszero(Symbolics.expand(s5 - s5_reference))
    end
  end
end


@testset "nonlinear oscillator" begin
  @testset "SymEngine.jl" begin
    dt = SymEngine.symbols("dt")
    u = SymEngine.symbols("u_1, u_2")
    f = [-u[2], u[1]] / (u[1]^2 + u[2]^2)

    # explicit midpoint method
    A = @SArray [0 0; 1//2 0]
    b = @SArray [0, 1//1]
    c = @SArray [0, 1//2];

    # tested with Mathematica using
    # ```
    # ClearAll["Global`*"];
    #
    # solution[t_,
    #   u0_] := {u0[[1]]*Cos[t/(u0[[1]]^2 + u0[[2]]^2)] -
    #   u0[[2]]*Sin[t/(u0[[1]]^2 + u0[[2]]^2)],
    #   u0[[2]]*Cos[t/(u0[[1]]^2 + u0[[2]]^2)] +
    #   u0[[1]]*Sin[t/(u0[[1]]^2 + u0[[2]]^2)]}
    #
    # factorF2[u_] := 1/(12*(u[[1]]^2 + u[[2]]^2)^2)
    # (*factorF2[u_]:=ff2*)
    # factorF4[u_] := 1/(20*(u[[1]]^2 + u[[2]]^2)^4)
    # (*factorF4[u_]:=ff4*)
    # factorF6[u_] := 127/(2016*(u[[1]]^2 + u[[2]]^2)^6)
    # (*factorF6[u_]:=ff6*)
    # factorF8[u_] := 1513/(12960*(u[[1]]^2 + u[[2]]^2)^8)
    # (*factorF8[u_]:=ff8*)
    #
    # factorU3[u_] := 0
    # (*factorU3[u_]:=fu3*)
    # factorU5[u_] := 1/(48*(u[[1]]^2 + u[[2]]^2)^6)
    # (*factorU5[u_]:=fu5*)
    # factorU7[u_] := 31/(640*(u[[1]]^2 + u[[2]]^2)^8)
    # (*factorU7[u_]:=fu7*)
    # factorU9[u_] := 8969/(80640*(u[[1]]^2 + u[[2]]^2)^10)
    # (*factorU9[u_]:=fu9*)
    #
    # newModifiedF[
    #   u_] :=
    # (1 + factorF2[u]*h^2 + factorF4[u]*h^4 + factorF6[u]*h^6)*
    #   originalF[
    #     u] + (factorU3[u]*h^3 + factorU5[u]*h^5 + factorU7[u]*h^7)*u
    #
    # u0 = {u01, u02};
    # f = newModifiedF;
    # y2 = u0 + h/2*f[u0];
    # unew = u0 + h*f[y2];
    # difference = Simplify@Series[unew - solution[h, u0], {h, 0, 10}];
    #
    # Simplify[difference[[1]]]
    # Simplify[difference[[2]]]
    # ```
    series = modifying_integrator(f, u, dt, A, b, c, 10)
    terms = SymEngine.subs.(series, (Dict(u[1] => 1//1, u[2] => 0//1), ))

    @test isequal(terms[1], 1//48 * dt^5 + 31//640 * dt^7 + 8969//80640 * dt^9)
    @test isequal(terms[2], 1 + 1//12 * dt^2 + 1//20 * dt^4 + 127//2016 * dt^6 + 1513//12960 * dt^8)
  end
end

end # @testset "BSeries"
