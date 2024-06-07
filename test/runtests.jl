using Test
using BSeries

using BSeries.Latexify: latexify

using LinearAlgebra: I
using StaticArrays: @SArray, @SMatrix, @SVector

using SymEngine: SymEngine
using SymPy: SymPy
using Symbolics: Symbolics

using Aqua: Aqua

@testset "BSeries" begin
    @testset "lazy representation of exact ODE solution" begin
        exact = ExactSolution{Rational{Int}}()
        terms = collect(Iterators.take(exact, 5))
        @test terms == [1 // 1, 1 // 1, 1 // 2, 1 // 6, 1 // 3]
        @test exact == ExactSolution(exact)
    end

    @testset "non-conflicting exports" begin
        # classical RK4
        A = [0 0 0 0
             1//2 0 0 0
             0 1//2 0 0
             0 0 1 0]
        b = [1 // 6, 1 // 3, 1 // 3, 1 // 6]
        rk = @inferred RungeKuttaMethod(A, b)
        t = @inferred rootedtree([1, 2])
        @test_nowarn @inferred derivative_weight(t, rk)
    end

    @testset "latexify" begin
        # explicit midpoint method
        A = @SArray [0 0; 1//2 0]
        b = @SArray [0, 1 // 1]
        c = @SArray [0, 1 // 2]

        series_integrator = bseries(A, b, c, 2)
        @test_nowarn latexify(series_integrator)
        @test_nowarn latexify(series_integrator, cdot = false)
        @test_nowarn latexify(series_integrator, dt = SymEngine.symbols("h"))
        @test_nowarn latexify(series_integrator - series_integrator)

        @testset "SymEngine.jl" begin
            @testset "Symbolic coefficients" begin
                α = SymEngine.symbols("α")
                A = [0 0; 1/(2 * α) 0]
                b = [1 - α, α]
                c = [0, 1 / (2 * α)]
                series_integrator = @inferred bseries(A, b, c, 3)
                @test_nowarn latexify(series_integrator)
            end

            @testset "Divide by h" begin
                A = [0 0; 1//2 0]
                b = [0, 1]
                c = [0, 1 // 2]
                series_integrator = @inferred bseries(A, b, c, 1)
                @test_nowarn latexify(series_integrator)

                h = SymEngine.symbols("h")
                coefficients = @inferred modified_equation(series_integrator)
                val1 = @test_nowarn latexify(coefficients, reduce_order_by = 1,
                                             cdot = false)
                val2 = @test_nowarn latexify(coefficients / h, cdot = false)
                @test val1 == val2
            end
        end

        @testset "SymPy.jl" begin
            @testset "Symbolic coefficients" begin
                α = SymPy.symbols("α", real = true)
                A = [0 0; 1/(2 * α) 0]
                b = [1 - α, α]
                c = [0, 1 / (2 * α)]
                series_integrator = bseries(A, b, c, 3)

                # Call this once to avoid failing tests with `@test_nowarn`
                # caused by deprecation warnings from PyCall.jl. See also
                # https://github.com/JuliaPy/PyCall.jl/pull/1042
                deepcopy(α)

                @test_nowarn latexify(series_integrator)
            end

            @testset "Divide by h" begin
                A = [0 0; 1//2 0]
                b = [0, 1]
                c = [0, 1 // 2]
                series_integrator = @inferred bseries(A, b, c, 1)
                @test_nowarn latexify(series_integrator)

                h = SymPy.symbols("h", real = true)
                coefficients = @inferred modified_equation(series_integrator)
                val1 = @test_nowarn latexify(coefficients, reduce_order_by = 1,
                                             cdot = false)
                val2 = @test_nowarn latexify(coefficients / h, cdot = false)
                @test val1 == val2
            end
        end

        @testset "Symbolics.jl" begin
            @testset "Symbolic coefficients" begin
                Symbolics.@variables α
                A = [0 0; 1/(2 * α) 0]
                b = [1 - α, α]
                c = [0, 1 / (2 * α)]
                series_integrator = @inferred bseries(A, b, c, 3)
                # Do not test warnings due to deprecation warnings from Symbolics
                # see https://github.com/ranocha/BSeries.jl/pull/210
                # @test_nowarn latexify(series_integrator)
                latexify(series_integrator)
            end

            @testset "Divide by h" begin
                A = [0 0; 1//2 0]
                b = [0, 1]
                c = [0, 1 // 2]
                series_integrator = @inferred bseries(A, b, c, 1)
                @test_nowarn latexify(series_integrator)

                Symbolics.@variables h
                coefficients = @inferred modified_equation(series_integrator)
                val1 = @test_nowarn latexify(coefficients, reduce_order_by = 1,
                                             cdot = false)
                val2 = @test_nowarn latexify(coefficients / h, cdot = false)
                @test val1 == val2
            end
        end
    end # @testset "latexify"

    @testset "AbstractDict interface" begin
        # explicit midpoint method
        A = @SArray [0 0; 1//2 0]
        b = @SArray [0, 1 // 1]
        c = @SArray [0, 1 // 2]
        series_integrator = bseries(A, b, c, 2)

        # These are mostly simple smoke tests
        @test_nowarn begin
            @test isempty(empty(series_integrator,
                                RootedTrees.RootedTree{Int32, Vector{Int32}},
                                Rational{Int32}))
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
    end # @testset "AbstractDict interface"

    @testset "vector space interface" begin
        # explicit midpoint method
        A = @SArray [0 0; 1//2 0]
        b = @SArray [0, 1 // 1]
        c = @SArray [0, 1 // 2]
        series1 = @inferred bseries(A, b, c, 2)
        series2 = @inferred bseries(float.(A), b, c, 3)
        exact = ExactSolution{Rational{Int128}}()

        diff = @inferred series1 - series2
        @test mapreduce(iszero, &, values(diff))

        @test @inferred(series1-series2) == @inferred(-(series2 - series1))

        diff = @inferred series1 - exact
        @test mapreduce(iszero, &, values(diff))

        diff = @inferred exact - series1
        @test mapreduce(iszero, &, values(diff))

        @test @inferred(series1+series2) == @inferred(2*series1)
        @test @inferred(series1+exact) == @inferred(2*series1)
        @test @inferred(exact+series1) == @inferred(2*series1)

        half1 = @inferred 0.5 * series1
        half2 = @inferred 2 \ series1
        half3 = @inferred series1 * 0.5
        half4 = @inferred series1 / 2
        @test half1 == half2
        @test half1 == half2
        @test half1 == half3
        @test half1 == half4

        @test @inferred(+series2) == series2
        @test @inferred(-series2) == -1 * series2

        diff = @inferred(series2+(-series2))
        @test mapreduce(iszero, &, values(diff))
    end # @testset "vector space interface"

    @testset "order_of_accuracy" begin
        @testset "RK4, rational coefficients" begin
            # classical fourth-order Runge-Kutta method
            A = [0 0 0 0;
                 1//2 0 0 0;
                 0 1//2 0 0;
                 0 0 1 0]
            b = [1 // 6, 1 // 3, 1 // 3, 1 // 6]
            rk = RungeKuttaMethod(A, b)
            series = bseries(rk, 5)
            @test @inferred(order_of_accuracy(series)) == 4

            series = bseries(rk, 3)
            @test @inferred(order_of_accuracy(series)) == 3
        end

        @testset "RK4, floating point coefficients" begin
            # classical fourth-order Runge-Kutta method
            A = [0 0 0 0;
                 1/2 0 0 0;
                 0 1/2 0 0;
                 0 0 1 0]
            b = [1 / 6, 1 / 3, 1 / 3, 1 / 6]
            rk = RungeKuttaMethod(A, b)
            series = bseries(rk, 5)
            @test @inferred(order_of_accuracy(series)) == 4
            @test @inferred(order_of_accuracy(series; atol = 10 * eps())) == 4

            series = bseries(rk, 3)
            @test @inferred(order_of_accuracy(series)) == 3
        end
    end

    @testset "average Vector field (AVF) method" begin
        @testset "Rational coefficients" begin
            series = @inferred(bseries(AverageVectorFieldMethod(), 6))
            @test @inferred(order_of_accuracy(series)) == 2
            @test is_energy_preserving(series)
        end

        @testset "Floating point coefficients" begin
            series = @inferred(bseries(AverageVectorFieldMethod(Float64), 6))
            @test @inferred(order_of_accuracy(series)) == 2
            @test is_energy_preserving(series)
        end
    end

    @testset "substitute" begin
        Symbolics.@variables a1 a2 a31 a32
        a = OrderedDict{RootedTrees.RootedTree{Int, Vector{Int}}, Symbolics.Num}(rootedtree(Int[]) => 1,
                                                                                 rootedtree([1]) => a1,
                                                                                 rootedtree([1, 2]) => a2,
                                                                                 rootedtree([1, 2, 2]) => a31,
                                                                                 rootedtree([1, 2, 3]) => a32)

        Symbolics.@variables b1 b2 b31 b32
        b = OrderedDict{RootedTrees.RootedTree{Int, Vector{Int}}, Symbolics.Num}(rootedtree(Int[]) => 0,
                                                                                 rootedtree([1]) => b1,
                                                                                 rootedtree([1, 2]) => b2,
                                                                                 rootedtree([1, 2, 2]) => b31,
                                                                                 rootedtree([1, 2, 3]) => b32)

        # See equation (6) of
        # - Philippe Chartier, Ernst Hairer and Gilles Vilmart (2007)
        #   Numerical integrators based on modified differential equations
        #   [DOI: 10.1090/S0025-5718-07-01967-9](https://doi.org/10.1090/S0025-5718-07-01967-9)

        b_a = @inferred substitute(b, a)
        t = rootedtree(Int[])
        @test isequal(b_a[t], a[t])

        t = rootedtree([1])
        coef = a1 * b1
        @inferred substitute(b, a, t)
        @test isequal(substitute(b, a, t), coef)
        @test isequal(b_a[t], coef)

        t = rootedtree([1, 2])
        coef = a1 * b2 + a2 * b1^2
        @inferred substitute(b, a, t)
        @test isequal(substitute(b, a, t), coef)
        @test isequal(b_a[t], coef)

        t = rootedtree([1, 2, 2])
        coef = a1 * b31 + 2 * a2 * b1 * b2 + a31 * b1^3
        @inferred substitute(b, a, t)
        @test isequal(substitute(b, a, t), coef)
        @test isequal(b_a[t], coef)

        t = rootedtree([1, 2, 3])
        coef = a1 * b32 + 2 * a2 * b1 * b2 + a32 * b1^3
        @inferred substitute(b, a, t)
        @test isequal(substitute(b, a, t), coef)
        @test isequal(b_a[t], coef)
    end # @testset "substitute"

    @testset "compose" begin
        Symbolics.@variables a0 a1 a2 a31 a32
        a = OrderedDict{RootedTrees.RootedTree{Int, Vector{Int}}, Symbolics.Num}(rootedtree(Int[]) => a0,
                                                                                 rootedtree([1]) => a1,
                                                                                 rootedtree([1, 2]) => a2,
                                                                                 rootedtree([1, 2, 2]) => a31,
                                                                                 rootedtree([1, 2, 3]) => a32)

        Symbolics.@variables b1 b2 b31 b32
        b = OrderedDict{RootedTrees.RootedTree{Int, Vector{Int}}, Symbolics.Num}(rootedtree(Int[]) => 0,
                                                                                 rootedtree([1]) => b1,
                                                                                 rootedtree([1, 2]) => b2,
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

        @testset "Composing RK4 with itself" begin
            # classical fourth-order Runge-Kutta method
            A = @SArray [0 0 0 0;
                         1//2 0 0 0;
                         0 1//2 0 0;
                         0 0 1 0]
            b = @SArray [1 // 6, 1 // 3, 1 // 3, 1 // 6]
            c = @SArray [0, 1 // 2, 1 // 2, 1]

            series_rk4 = bseries(A, b, c, 8)
            series_rk4_composed = @inferred compose(series_rk4, series_rk4,
                                                    normalize_stepsize = true)

            # Butcher coefficients of the composition (with normalized step size)
            A = @SArray [0 0 0 0 0 0 0 0;
                         1//4 0 0 0 0 0 0 0;
                         0 1//4 0 0 0 0 0 0;
                         0 0 1//2 0 0 0 0 0;
                         1//12 1//6 1//6 1//12 0 0 0 0;
                         1//12 1//6 1//6 1//12 1//4 0 0 0;
                         1//12 1//6 1//6 1//12 0 1//4 0 0;
                         1//12 1//6 1//6 1//12 0 0 1//2 0]
            b = @SArray [1 // 12, 1 // 6, 1 // 6, 1 // 12, 1 // 12, 1 // 6, 1 // 6, 1 // 12]
            c = @SArray [0, 1 // 4, 1 // 4, 1 // 2, 1 // 2, 3 // 4, 3 // 4, 1 // 1]

            series_2rk4 = bseries(A, b, c, order(series_rk4))

            @test series_rk4_composed == series_2rk4

            # Ensure that the normalization uses the correct factor
            unnormalized_series = @inferred compose(series_rk4, series_rk4,
                                                    normalize_stepsize = false)
            for t in keys(unnormalized_series)
                @test unnormalized_series[t] == series_rk4_composed[t] * 2^order(t)
            end
        end

        @testset "Butcher's effective order 5 method" begin
            # Butcher, J. C.
            # "The effective order of Runge-Kutta methods."
            # In Conference on the numerical solution of differential equations,
            # pp. 133-139. Springer, Berlin, Heidelberg, 1969.
            # https://doi.org/10.1007/BFb0060019
            A = [0 0 0 0 0;
                 1//5 0 0 0 0;
                 0 2//5 0 0 0;
                 3//16 0 5//16 0 0;
                 1//4 0 -5//4 2 0]
            b = [1 // 6, 0, 0, 2 // 3, 1 // 6]
            rk_a = RungeKuttaMethod(A, b)
            series_a = @inferred bseries(rk_a, 5)
            # this is the main method

            A = [0 0 0 0 0;
                 1//5 0 0 0 0;
                 0 2//5 0 0 0;
                 75//64 -9//4 117//64 0 0;
                 -37//36 7//3 -3//4 4//9 0]
            b = [19 // 144, 0, 25 // 48, 2 // 9, 1 // 8]
            rk_b = RungeKuttaMethod(A, b)
            series_b = @inferred bseries(rk_b, 5)
            # this is the starting method

            A = [0 0 0 0 0;
                 1//5 0 0 0 0;
                 0 2//5 0 0 0;
                 161//192 -19//12 287//192 0 0;
                 -27//28 19//7 -291//196 36//49 0]
            b = [7 // 48, 0, 475 // 1008, 2 // 7, 7 // 72]
            rk_c = RungeKuttaMethod(A, b)
            series_c = @inferred bseries(rk_c, 5)
            # this is the finishing method

            series = @inferred compose(series_b, series_a, series_c,
                                       normalize_stepsize = true)
            @test series == ExactSolution(series)
        end
    end # @testset "compose"

    @testset "modified_equation" begin
        t1 = rootedtree([1])
        t2 = rootedtree([1, 2])
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
            b = @SArray [1 / 2, 1 / 2]
            c = @SArray [0.0, 1.0]
            order = 5
            series = modified_equation(A, b, c, order)

            # tested with the Python package BSeries
            @test series[t1]≈1.0 atol=10 * eps()
            @test series[t2]≈0 atol=10 * eps()
            @test series[t31]≈0.166666666666667 atol=10 * eps()
            @test series[t32]≈-0.166666666666667 atol=10 * eps()
            @test series[t41]≈-1.11022302462516e-16 atol=10 * eps()
            @test series[t42]≈-0.125000000000000 atol=10 * eps()
            @test series[t43]≈-2.77555756156289e-17 atol=10 * eps()
            @test series[t44]≈0.125000000000000 atol=10 * eps()
            @test series[t51]≈-0.0333333333333332 atol=10 * eps()
            @test series[t52]≈-0.0583333333333333 atol=10 * eps()
            @test series[t55]≈0.0750000000000000 atol=10 * eps()
            @test series[t53]≈0.0583333333333333 atol=10 * eps()
            @test series[t56]≈0.0333333333333334 atol=10 * eps()
            @test series[t54]≈0.0500000000000000 atol=10 * eps()
            @test series[t57]≈0.0583333333333333 atol=10 * eps()
            @test series[t58]≈-0.0583333333333333 atol=10 * eps()
            @test series[t59]≈-0.0500000000000000 atol=10 * eps()
        end

        @testset "SSPRK(3, 3)" begin
            A = @SArray [0.0 0.0 0.0;
                         1.0 0.0 0.0;
                         1/4 1/4 0.0]
            b = @SArray [1 / 6, 1 / 6, 2 / 3]
            c = @SArray [0.0, 1.0, 1 / 2]
            order = 5
            series = modified_equation(A, b, c, order)

            # tested with the Python package BSeries
            @test series[t1]≈1.0 atol=10 * eps()
            @test series[t2]≈0 atol=10 * eps()
            @test series[t31]≈0 atol=10 * eps()
            @test series[t32]≈0 atol=10 * eps()
            @test series[t41]≈0 atol=10 * eps()
            @test series[t42]≈-0.0416666666666667 atol=10 * eps()
            @test series[t43]≈0.0833333333333333 atol=10 * eps()
            @test series[t44]≈-0.0416666666666667 atol=10 * eps()
            @test series[t51]≈0.00833333333333330 atol=10 * eps()
            @test series[t52]≈-0.0166666666666667 atol=10 * eps()
            @test series[t55]≈0.0333333333333333 atol=10 * eps()
            @test series[t53]≈0.0166666666666667 atol=10 * eps()
            @test series[t56]≈-0.00833333333333333 atol=10 * eps()
            @test series[t54]≈0.00833333333333333 atol=10 * eps()
            @test series[t57]≈-0.0250000000000000 atol=10 * eps()
            @test series[t58]≈-0.0166666666666667 atol=10 * eps()
            @test series[t59]≈0.0333333333333333 atol=10 * eps()
        end

        @testset "Explicit midpoint method (Runge's method)" begin
            A = @SArray [0 0; 1//2 0]
            b = @SArray [0, 1 // 1]
            c = @SArray [0, 1 // 2]
            order = 5
            series = modified_equation(A, b, c, order)

            # tested with the Python package BSeries
            @test series[t1] == 1
            @test series[t2] == 0
            @test series[t31] == -1 // 12
            @test series[t32] == -1 // 6
            @test series[t41] == 0
            @test series[t42] == 0
            @test series[t43] == 1 // 8
            @test series[t44] == 1 // 8
            @test series[t51] == 7 // 240
            @test series[t52] == 1 // 40
            @test series[t55] == 1 // 30
            @test series[t53] == 3 // 80
            @test series[t56] == -7 // 240
            @test series[t54] == 7 // 240
            @test series[t57] == -1 // 40
            @test series[t58] == -19 // 240
            @test series[t59] == -1 // 20
        end

        @testset "Average vector field method" begin
            # Example 1 of
            # Elena Celledoni, Robert I. McLachlan, Brynjulf Owren, and G. R. W. Quispel (2010)
            # Energy-preserving integrators and the structure of B-series.
            # Foundations of Computational Mathematics
            # [DOI: 10.1007/s10208-010-9073-1](https://doi.org/10.1007/s10208-010-9073-1)
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

            coefficients = @inferred modified_equation(series)
            let t = rootedtree(Int[])
                @test coefficients[t] / symmetry(t) == 0
            end
            let t = rootedtree([1])
                @test coefficients[t] / symmetry(t) == 1
            end
            let t = rootedtree([1, 2])
                @test coefficients[t] / symmetry(t) == 0
            end
            let t = rootedtree([1, 2, 3])
                @test coefficients[t] / symmetry(t) == 1 // 12
            end
            let t = rootedtree([1, 2, 2])
                @test coefficients[t] / symmetry(t) == 0
            end
            let t = rootedtree([1, 2, 3, 4])
                @test coefficients[t] / symmetry(t) == 0
            end
            let t = rootedtree([1, 2, 3, 3])
                @test coefficients[t] / symmetry(t) == 0
            end
            let t = rootedtree([1, 2, 3, 2])
                @test coefficients[t] / symmetry(t) == 0
            end
            let t = rootedtree([1, 2, 2, 2])
                @test coefficients[t] / symmetry(t) == 0
            end
            let t = rootedtree([1, 2, 3, 4, 5])
                @test coefficients[t] / symmetry(t) == 9 // 720
            end
            let t = rootedtree([1, 2, 3, 4, 4])
                @test coefficients[t] / symmetry(t) == 4 // 720
            end
            let t = rootedtree([1, 2, 3, 4, 3])
                @test coefficients[t] / symmetry(t) == 2 // 720
            end
            let t = rootedtree([1, 2, 3, 4, 2])
                @test coefficients[t] / symmetry(t) == -4 // 720
            end
            let t = rootedtree([1, 2, 3, 3, 3])
                @test coefficients[t] / symmetry(t) == -1 // 720
            end
            let t = rootedtree([1, 2, 3, 3, 2])
                @test coefficients[t] / symmetry(t) == -4 // 720
            end
            let t = rootedtree([1, 2, 3, 2, 3])
                @test coefficients[t] / symmetry(t) == 2 // 720
            end
            let t = rootedtree([1, 2, 3, 2, 2])
                @test coefficients[t] / symmetry(t) == -1 // 720
            end
            let t = rootedtree([1, 2, 2, 2, 2])
                @test coefficients[t] / symmetry(t) == 0
            end
        end
    end # @testset "modified_equation"

    @testset "elementary differentials" begin
        @testset "Bicolored trees" begin
            @testset "Lotka-Volterra" begin
                # Verified with Mathematica
                # u = {q, p};
                # f1 = {q * (p - 1), 0};
                # f2 = {0, p * (2 - q)};
                #
                # f1p = Simplify@D[f1, {u}];
                # f2p = Simplify@D[f2, {u}];
                #
                # f1pp = Simplify@D[f1, {u, 2}];
                # f2pp = Simplify@D[f2, {u, 2}];
                #
                # (* Trees of order 2 *)
                # Map[InputForm, f1p.f1]
                # Map[InputForm, f2p.f1]
                # Map[InputForm, f1p.f2]
                # Map[InputForm, f2p.f2]
                #
                # (* Tall trees of order 3 *)
                # Map[InputForm, f1p.(f1p.f1)]
                # Map[InputForm, f2p.(f1p.f1)]
                # Map[InputForm, f1p.(f2p.f1)]
                # Map[InputForm, f2p.(f2p.f1)]
                # Map[InputForm, f1p.(f1p.f2)]
                # Map[InputForm, f2p.(f1p.f2)]
                # Map[InputForm, f1p.(f2p.f2)]
                # Map[InputForm, f2p.(f2p.f2)]
                #
                # (* Bushy trees of order 3 *)
                # Map[InputForm, f1pp.f1.f1]
                # Map[InputForm, f2pp.f1.f1]
                # Map[InputForm, f1pp.f2.f1]
                # Map[InputForm, f2pp.f2.f1]
                # (* The other choice of colors of the leaves is redundant *)
                # Map[InputForm, f1pp.f2.f2]
                # Map[InputForm, f2pp.f2.f2]
                @testset "SymEngine" begin
                    q, p = u = SymEngine.symbols("q, p")
                    f = ([q * (p - 1), 0], [0, p * (2 - q)])
                    differentials = elementary_differentials(f, u, 3)
                    @test differentials[rootedtree(Int64[], Bool[])] == [1, 1]
                    @test differentials[rootedtree([1], Bool[0])] == [q * (p - 1), 0]
                    @test differentials[rootedtree([1], Bool[1])] == [0, p * (2 - q)]
                    @test differentials[rootedtree([1, 2], Bool[0, 0])] ==
                          [(-1 + p)^2 * q, 0]
                    @test differentials[rootedtree([1, 2], Bool[1, 0])] ==
                          [0, -((-1 + p) * p * q)]
                    @test differentials[rootedtree([1, 2], Bool[0, 1])] ==
                          [p * (2 - q) * q, 0]
                    @test differentials[rootedtree([1, 2], Bool[1, 1])] ==
                          [0, p * (2 - q)^2]
                    @test differentials[rootedtree([1, 2, 3], Bool[0, 0, 0])] ==
                          [(-1 + p)^3 * q, 0]
                    @test differentials[rootedtree([1, 2, 3], Bool[1, 0, 0])] ==
                          [0, -((-1 + p)^2 * p * q)]
                    @test differentials[rootedtree([1, 2, 3], Bool[0, 1, 0])] ==
                          [-((-1 + p) * p * q^2), 0]
                    @test differentials[rootedtree([1, 2, 3], Bool[1, 1, 0])] ==
                          [0, -((-1 + p) * p * (2 - q) * q)]
                    @test differentials[rootedtree([1, 2, 3], Bool[0, 0, 1])] ==
                          [(-1 + p) * p * (2 - q) * q, 0]
                    @test differentials[rootedtree([1, 2, 3], Bool[1, 0, 1])] ==
                          [0, -(p^2 * (2 - q) * q)]
                    @test differentials[rootedtree([1, 2, 3], Bool[0, 1, 1])] ==
                          [p * (2 - q)^2 * q, 0]
                    @test differentials[rootedtree([1, 2, 3], Bool[1, 1, 1])] ==
                          [0, p * (2 - q)^3]
                    @test differentials[rootedtree([1, 2, 2], Bool[0, 0, 0])] == [0, 0]
                    @test differentials[rootedtree([1, 2, 2], Bool[1, 0, 0])] == [0, 0]
                    @test differentials[rootedtree([1, 2, 2], Bool[0, 1, 0])] ==
                          [(-1 + p) * p * (2 - q) * q, 0]
                    @test differentials[rootedtree([1, 2, 2], Bool[1, 1, 0])] ==
                          [0, -((-1 + p) * p * (2 - q) * q)]
                    @test differentials[rootedtree([1, 2, 2], Bool[0, 1, 1])] == [0, 0]
                    @test differentials[rootedtree([1, 2, 2], Bool[1, 1, 1])] == [0, 0]
                end
            end
        end
    end # @testset "elementary differentials"

    @testset "modified_equation with elementary differentials" begin
        @testset "SymEngine.jl" begin
            @testset "Explicit Euler" begin
                # Lotka-Volterra model
                dt = SymEngine.symbols("dt")
                p, q = u = SymEngine.symbols("p, q")
                f = [p * (2 - q), q * (p - 1)]

                # Explicit Euler method
                A = @SMatrix [0 // 1;;]
                b = @SArray [1 // 1]
                c = @SArray [0 // 1]

                # tested with the Python package BSeries
                series2 = modified_equation(f, u, dt, A, b, c, 2)
                series2_reference = [
                    -dt * (-p * q * (p - 1) + p * (2 - q)^2) / 2 + p * (2 - q),
                    -dt * (p * q * (2 - q) + q * (p - 1)^2) / 2 + q * (p - 1),
                ]
                @test mapreduce(isequal, &, series2, series2_reference)

                # tested with the Python package BSeries
                series3 = modified_equation(f, u, dt, A, b, c, 3)
                series3_reference = [
                    -dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                    dt^2 *
                    (-p^2 * q * (2 - q) - p * q * (2 - q) * (p - 1) - p * q * (p - 1)^2 +
                     p * (2 - q)^3) / 3 - dt * (-p * q * (p - 1) + p * (2 - q)^2) / 2 +
                    p * (2 - q),
                    dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                    dt^2 *
                    (-p * q^2 * (p - 1) + p * q * (2 - q)^2 + p * q * (2 - q) * (p - 1) +
                     q * (p - 1)^3) / 3 - dt * (p * q * (2 - q) + q * (p - 1)^2) / 2 +
                    q * (p - 1),
                ]
                @test mapreduce(iszero ∘ SymEngine.expand, &, series3 - series3_reference)

                # tested with the Python package BSeries
                series4 = modified_equation(f, u, dt, A, b, c, 4)
                series4_reference = [
                    -dt^3 * (-2 * p^2 * q * (2 - q) * (p - 1) +
                     2 * p * q * (2 - q) * (p - 1) * (q - 2)) / 12 -
                    dt^3 * (-p^2 * q * (2 - q)^2 + p * q^2 * (p - 1)^2 +
                     p * q * (1 - p) * (2 - q) * (p - 1) +
                     p * q * (2 - q) * (p - 1) * (q - 2)) / 12 -
                    dt^3 * (p^2 * q^2 * (p - 1) - 2 * p^2 * q * (2 - q)^2 -
                     p^2 * q * (2 - q) * (p - 1) - p * q * (2 - q)^2 * (p - 1) -
                     p * q * (2 - q) * (p - 1)^2 - p * q * (p - 1)^3 + p * (2 - q)^4) / 4 -
                    dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                    dt^2 *
                    (-p^2 * q * (2 - q) - p * q * (2 - q) * (p - 1) - p * q * (p - 1)^2 +
                     p * (2 - q)^3) / 3 - dt * (-p * q * (p - 1) + p * (2 - q)^2) / 2 +
                    p * (2 - q),
                    -dt^3 *
                    (-2 * p * q^2 * (2 - q) * (p - 1) + 2 * p * q * (2 - q) * (p - 1)^2) /
                    12 -
                    dt^3 * (p^2 * q * (2 - q)^2 - p * q^2 * (p - 1)^2 +
                     p * q * (2 - q)^2 * (p - 1) + p * q * (2 - q) * (p - 1)^2) / 12 -
                    dt^3 * (-p^2 * q^2 * (2 - q) - p * q^2 * (2 - q) * (p - 1) -
                     2 * p * q^2 * (p - 1)^2 + p * q * (2 - q)^3 +
                     p * q * (2 - q)^2 * (p - 1) + p * q * (2 - q) * (p - 1)^2 +
                     q * (p - 1)^4) / 4 + dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                    dt^2 *
                    (-p * q^2 * (p - 1) + p * q * (2 - q)^2 + p * q * (2 - q) * (p - 1) +
                     q * (p - 1)^3) / 3 - dt * (p * q * (2 - q) + q * (p - 1)^2) / 2 +
                    q * (p - 1),
                ]
                @test mapreduce(iszero ∘ SymEngine.expand, &, series4 - series4_reference)
            end

            @testset "IMEX Euler (partitioned Euler, symplectic Euler)" begin
                # Lotka-Volterra model
                dt = SymEngine.symbols("dt")
                q, p = u = SymEngine.symbols("q, p")
                f = ([0, p * (2 - q)], [q * (p - 1), 0])

                # Symplectic Euler method
                ex_euler = RungeKuttaMethod(@SMatrix([0]), @SVector [1])
                im_euler = RungeKuttaMethod(@SMatrix([1]), @SVector [1])
                ark = AdditiveRungeKuttaMethod([im_euler, ex_euler])

                # Hairer, Lubich, Wanner (2006) Geometric numerical integration
                # Example IX.1.3(b)
                series_integrator = bseries(ark, 2)
                series = modified_equation(f, u, dt, series_integrator)
                series_reference = [
                    q * (p - 1) - dt / 2 * q * (p^2 + p * q - 4 * p + 1),
                    -p * (q - 2) + dt / 2 * p * (q^2 + p * q - 5 * q + 4),
                ]
                @test mapreduce(iszero ∘ SymEngine.expand, &, series - series_reference)
            end
        end

        @testset "SymPy.jl" begin
            # Lotka-Volterra model
            dt = SymPy.symbols("dt")
            p, q = u = SymPy.symbols("p, q")
            f = [p * (2 - q), q * (p - 1)]

            # Explicit Euler method
            A = @SMatrix [0 // 1;;]
            b = @SArray [1 // 1]
            c = @SArray [0 // 1]

            # tested with the Python package BSeries
            series2 = modified_equation(f, u, dt, A, b, c, 2)
            series2_reference = [
                -dt * (-p * q * (p - 1) + p * (2 - q)^2) / 2 + p * (2 - q),
                -dt * (p * q * (2 - q) + q * (p - 1)^2) / 2 + q * (p - 1),
            ]
            @test mapreduce(isequal, &, series2, series2_reference)

            # tested with the Python package BSeries
            series3 = modified_equation(f, u, dt, A, b, c, 3)
            series3_reference = [
                -dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                dt^2 *
                (-p^2 * q * (2 - q) - p * q * (2 - q) * (p - 1) - p * q * (p - 1)^2 +
                 p * (2 - q)^3) / 3 - dt * (-p * q * (p - 1) + p * (2 - q)^2) / 2 +
                p * (2 - q),
                dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                dt^2 *
                (-p * q^2 * (p - 1) + p * q * (2 - q)^2 + p * q * (2 - q) * (p - 1) +
                 q * (p - 1)^3) / 3 - dt * (p * q * (2 - q) + q * (p - 1)^2) / 2 +
                q * (p - 1),
            ]
            @test mapreduce(iszero ∘ SymPy.expand, &, series3 - series3_reference)

            # tested with the Python package BSeries
            series4 = modified_equation(f, u, dt, A, b, c, 4)
            series4_reference = [
                -dt^3 *
                (-2 * p^2 * q * (2 - q) * (p - 1) + 2 * p * q * (2 - q) * (p - 1) * (q - 2)) /
                12 -
                dt^3 * (-p^2 * q * (2 - q)^2 + p * q^2 * (p - 1)^2 +
                 p * q * (1 - p) * (2 - q) * (p - 1) + p * q * (2 - q) * (p - 1) * (q - 2)) /
                12 -
                dt^3 * (p^2 * q^2 * (p - 1) - 2 * p^2 * q * (2 - q)^2 -
                 p^2 * q * (2 - q) * (p - 1) - p * q * (2 - q)^2 * (p - 1) -
                 p * q * (2 - q) * (p - 1)^2 - p * q * (p - 1)^3 + p * (2 - q)^4) / 4 -
                dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                dt^2 *
                (-p^2 * q * (2 - q) - p * q * (2 - q) * (p - 1) - p * q * (p - 1)^2 +
                 p * (2 - q)^3) / 3 - dt * (-p * q * (p - 1) + p * (2 - q)^2) / 2 +
                p * (2 - q),
                -dt^3 *
                (-2 * p * q^2 * (2 - q) * (p - 1) + 2 * p * q * (2 - q) * (p - 1)^2) / 12 -
                dt^3 *
                (p^2 * q * (2 - q)^2 - p * q^2 * (p - 1)^2 + p * q * (2 - q)^2 * (p - 1) +
                 p * q * (2 - q) * (p - 1)^2) / 12 -
                dt^3 * (-p^2 * q^2 * (2 - q) - p * q^2 * (2 - q) * (p - 1) -
                 2 * p * q^2 * (p - 1)^2 + p * q * (2 - q)^3 + p * q * (2 - q)^2 * (p - 1) +
                 p * q * (2 - q) * (p - 1)^2 + q * (p - 1)^4) / 4 +
                dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                dt^2 *
                (-p * q^2 * (p - 1) + p * q * (2 - q)^2 + p * q * (2 - q) * (p - 1) +
                 q * (p - 1)^3) / 3 - dt * (p * q * (2 - q) + q * (p - 1)^2) / 2 +
                q * (p - 1),
            ]
            @test mapreduce(iszero ∘ SymPy.expand, &, series4 - series4_reference)
        end

        @testset "Symbolics.jl" begin
            # Lotka-Volterra model
            Symbolics.@variables dt
            u = Symbolics.@variables p q
            f = [p * (2 - q), q * (p - 1)]

            # Explicit Euler method
            A = @SMatrix [0 // 1;;]
            b = @SArray [1 // 1]
            c = @SArray [0 // 1]

            # tested with the Python package BSeries
            series2 = modified_equation(f, u, dt, A, b, c, 2)
            series2_reference = [
                -dt * (-p * q * (p - 1) + p * (2 - q)^2) / 2 + p * (2 - q),
                -dt * (p * q * (2 - q) + q * (p - 1)^2) / 2 + q * (p - 1),
            ]
            @test mapreduce(isequal, &, series2, series2_reference)

            # tested with the Python package BSeries
            series3 = modified_equation(f, u, dt, A, b, c, 3)
            series3_reference = [
                -dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                dt^2 *
                (-p^2 * q * (2 - q) - p * q * (2 - q) * (p - 1) - p * q * (p - 1)^2 +
                 p * (2 - q)^3) / 3 - dt * (-p * q * (p - 1) + p * (2 - q)^2) / 2 +
                p * (2 - q),
                dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                dt^2 *
                (-p * q^2 * (p - 1) + p * q * (2 - q)^2 + p * q * (2 - q) * (p - 1) +
                 q * (p - 1)^3) / 3 - dt * (p * q * (2 - q) + q * (p - 1)^2) / 2 +
                q * (p - 1),
            ]
            @test mapreduce(iszero ∘ Symbolics.expand, &, series3 - series3_reference)

            # tested with the Python package BSeries
            series4 = modified_equation(f, u, dt, A, b, c, 4)
            series4_reference = [
                -dt^3 *
                (-2 * p^2 * q * (2 - q) * (p - 1) + 2 * p * q * (2 - q) * (p - 1) * (q - 2)) /
                12 -
                dt^3 * (-p^2 * q * (2 - q)^2 + p * q^2 * (p - 1)^2 +
                 p * q * (1 - p) * (2 - q) * (p - 1) + p * q * (2 - q) * (p - 1) * (q - 2)) /
                12 -
                dt^3 * (p^2 * q^2 * (p - 1) - 2 * p^2 * q * (2 - q)^2 -
                 p^2 * q * (2 - q) * (p - 1) - p * q * (2 - q)^2 * (p - 1) -
                 p * q * (2 - q) * (p - 1)^2 - p * q * (p - 1)^3 + p * (2 - q)^4) / 4 -
                dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                dt^2 *
                (-p^2 * q * (2 - q) - p * q * (2 - q) * (p - 1) - p * q * (p - 1)^2 +
                 p * (2 - q)^3) / 3 - dt * (-p * q * (p - 1) + p * (2 - q)^2) / 2 +
                p * (2 - q),
                -dt^3 *
                (-2 * p * q^2 * (2 - q) * (p - 1) + 2 * p * q * (2 - q) * (p - 1)^2) / 12 -
                dt^3 *
                (p^2 * q * (2 - q)^2 - p * q^2 * (p - 1)^2 + p * q * (2 - q)^2 * (p - 1) +
                 p * q * (2 - q) * (p - 1)^2) / 12 -
                dt^3 * (-p^2 * q^2 * (2 - q) - p * q^2 * (2 - q) * (p - 1) -
                 2 * p * q^2 * (p - 1)^2 + p * q * (2 - q)^3 + p * q * (2 - q)^2 * (p - 1) +
                 p * q * (2 - q) * (p - 1)^2 + q * (p - 1)^4) / 4 +
                dt^2 * p * q * (2 - q) * (p - 1) / 6 +
                dt^2 *
                (-p * q^2 * (p - 1) + p * q * (2 - q)^2 + p * q * (2 - q) * (p - 1) +
                 q * (p - 1)^3) / 3 - dt * (p * q * (2 - q) + q * (p - 1)^2) / 2 +
                q * (p - 1),
            ]
            @test mapreduce(iszero ∘ Symbolics.expand, &, series4 - series4_reference)
        end
    end # @testset "modified_equation with elementary differentials"

    @testset "modifying_integrator" begin
        t1 = rootedtree([1])
        t2 = rootedtree([1, 2])
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
            b = @SArray [1 / 2, 1 / 2]
            c = @SArray [0.0, 1.0]
            order = 5
            series = modifying_integrator(A, b, c, order)

            # tested with the Python package BSeries
            @test series[t1]≈1 atol=10 * eps()
            @test series[t2]≈0 atol=10 * eps()
            @test series[t31]≈-0.166666666666667 atol=10 * eps()
            @test series[t32]≈1 / 6 atol=10 * eps()
            @test series[t41]≈8.32667268468867e-17 atol=10 * eps()
            @test series[t42]≈0.125000000000000 atol=10 * eps()
            @test series[t43]≈1.38777878078145e-17 atol=10 * eps()
            @test series[t44]≈-0.125000000000000 atol=10 * eps()
            @test series[t51]≈0.200000000000000 atol=10 * eps()
            @test series[t52]≈0.0583333333333333 atol=10 * eps()
            @test series[t55]≈0.00833333333333335 atol=10 * eps()
            @test series[t53]≈-0.0583333333333333 atol=10 * eps()
            @test series[t56]≈-0.200000000000000 atol=10 * eps()
            @test series[t54]≈-0.133333333333333 atol=10 * eps()
            @test series[t57]≈-0.0583333333333333 atol=10 * eps()
            @test series[t58]≈0.0583333333333333 atol=10 * eps()
            @test series[t59]≈0.133333333333333 atol=10 * eps()
        end

        @testset "SSPRK(3, 3)" begin
            A = @SArray [0.0 0.0 0.0;
                         1.0 0.0 0.0;
                         1/4 1/4 0.0]
            b = @SArray [1 / 6, 1 / 6, 2 / 3]
            c = @SArray [0.0, 1.0, 1 / 2]
            order = 5
            series = modifying_integrator(A, b, c, order)

            # tested with the Python package BSeries
            @test series[t1]≈1 atol=10 * eps()
            @test series[t2]≈0 atol=10 * eps()
            @test series[t31]≈0 atol=10 * eps()
            @test series[t32]≈0 atol=10 * eps()
            @test series[t41]≈0 atol=10 * eps()
            @test series[t42]≈0.0416666666666667 atol=10 * eps()
            @test series[t43]≈-0.0833333333333333 atol=10 * eps()
            @test series[t44]≈1 / 24 atol=10 * eps()
            @test series[t51]≈-0.00833333333333330 atol=10 * eps()
            @test series[t52]≈0.0166666666666667 atol=10 * eps()
            @test series[t55]≈-0.0333333333333333 atol=10 * eps()
            @test series[t53]≈-0.0166666666666667 atol=10 * eps()
            @test series[t56]≈0.00833333333333332 atol=10 * eps()
            @test series[t54]≈-0.00833333333333334 atol=10 * eps()
            @test series[t57]≈0.0250000000000000 atol=10 * eps()
            @test series[t58]≈1 / 60 atol=10 * eps()
            @test series[t59]≈-0.0333333333333333 atol=10 * eps()
        end

        @testset "Explicit midpoint method (Runge's method)" begin
            A = @SArray [0 0; 1//2 0]
            b = @SArray [0, 1 // 1]
            c = @SArray [0, 1 // 2]
            order = 5
            series = modifying_integrator(A, b, c, order)

            # tested with the Python package BSeries
            @test series[t1] == 1
            @test series[t2] == 0
            @test series[t31] == 1 // 12
            @test series[t32] == 1 // 6
            @test series[t41] == 0
            @test series[t42] == 0
            @test series[t43] == -1 // 8
            @test series[t44] == -1 // 8
            @test series[t51] == 1 // 80
            @test series[t52] == 1 // 60
            @test series[t55] == 7 // 240
            @test series[t53] == 1 // 240
            @test series[t56] == 9 // 80
            @test series[t54] == 1 // 80
            @test series[t57] == 13 // 120
            @test series[t58] == 13 // 80
            @test series[t59] == 2 // 15
        end
    end # @testset "modifying_integrator"

    @testset "modifying_integrator with elementary differentials" begin
        @testset "SymEngine.jl" begin
            dt = SymEngine.symbols("dt")
            α, β, γ = SymEngine.symbols("alpha, beta, gamma")
            u1, u2, u3 = u = SymEngine.symbols("u1 u2 u3")
            f = [α * u[2] * u[3], β * u[3] * u[1], γ * u[1] * u[2]]

            # Implicit midpoint method
            A = @SMatrix [1 // 2;;]
            b = @SArray [1 // 1]
            c = @SArray [1 // 2]

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
                s3 = SymEngine.subs(1 // 2 *
                                    SymEngine.diff(SymEngine.diff((series[i] - f[i]) / f[i],
                                                                  dt), dt),
                                    dt => 0)
                @test iszero(SymEngine.expand(s3 - s3_reference))
            end

            s5_reference = 6 // 5 * s3_reference^2 +
                           1 // 60 * α * β * γ *
                           (β * u1^2 * u3^2 + γ * u2^2 * u1^2 + α * u3^2 * u2^2)
            for i in eachindex(f)
                s5 = SymEngine.subs(1 // 24 *
                                    SymEngine.diff(SymEngine.diff(SymEngine.diff(SymEngine.diff((series[i] -
                                                                                                 f[i]) /
                                                                                                f[i],
                                                                                                dt),
                                                                                 dt), dt),
                                                   dt),
                                    dt => 0)
                @test iszero(SymEngine.expand(s5 - s5_reference))
            end
        end

        @testset "SymPy.jl" begin
            dt = SymPy.symbols("dt")
            α, β, γ = SymPy.symbols("alpha, beta, gamma")
            u1, u2, u3 = u = SymPy.symbols("u1 u2 u3")
            f = [α * u[2] * u[3], β * u[3] * u[1], γ * u[1] * u[2]]

            # Implicit midpoint method
            A = @SMatrix [1 // 2;;]
            b = @SArray [1 // 1]
            c = @SArray [1 // 2]

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
                s3 = SymPy.subs(1 // 2 *
                                SymPy.diff(SymPy.diff((series[i] - f[i]) / f[i], dt), dt),
                                dt => 0)
                @test iszero(SymPy.expand(s3 - s3_reference))
            end

            s5_reference = 6 // 5 * s3_reference^2 +
                           1 // 60 * α * β * γ *
                           (β * u1^2 * u3^2 + γ * u2^2 * u1^2 + α * u3^2 * u2^2)
            for i in eachindex(f)
                s5 = SymPy.subs(1 // 24 *
                                SymPy.diff(SymPy.diff(SymPy.diff(SymPy.diff((series[i] -
                                                                             f[i]) / f[i],
                                                                            dt), dt), dt),
                                           dt),
                                dt => 0)
                @test iszero(SymPy.expand(s5 - s5_reference))
            end
        end

        @testset "Symbolics.jl" begin
            Symbolics.@variables dt
            Symbolics.@variables α β γ
            u = Symbolics.@variables u1 u2 u3
            f = [α * u[2] * u[3], β * u[3] * u[1], γ * u[1] * u[2]]

            # Implicit midpoint method
            A = @SMatrix [1 // 2;;]
            b = @SArray [1 // 1]
            c = @SArray [1 // 2]

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
                s3 = Symbolics.simplify_fractions(Symbolics.substitute(Symbolics.expand_derivatives(1 //
                                                                                                    2 *
                                                                                                    d(d((series[i] -
                                                                                                         f[i]) /
                                                                                                        f[i]))),
                                                                       dt => 0))
                @test iszero(Symbolics.expand(s3 - s3_reference))
            end

            s5_reference = 6 // 5 * s3_reference^2 +
                           1 // 60 * α * β * γ *
                           (β * u1^2 * u3^2 + γ * u2^2 * u1^2 + α * u3^2 * u2^2)
            for i in eachindex(f)
                s5 = Symbolics.simplify_fractions(Symbolics.substitute(Symbolics.expand_derivatives(1 //
                                                                                                    24 *
                                                                                                    d(d(d(d((series[i] -
                                                                                                             f[i]) /
                                                                                                            f[i]))))),
                                                                       dt => 0))
                @test iszero(Symbolics.expand(s5 - s5_reference))
            end
        end
    end # @testset "modifying_integrator with elementary differentials"

    @testset "nonlinear oscillator" begin
        @testset "SymEngine.jl" begin
            dt = SymEngine.symbols("dt")
            u = SymEngine.symbols("u_1, u_2")
            f = [-u[2], u[1]] / (u[1]^2 + u[2]^2)

            # explicit midpoint method
            A = @SArray [0 0; 1//2 0]
            b = @SArray [0, 1 // 1]
            c = @SArray [0, 1 // 2]

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
            terms = SymEngine.subs.(series, (Dict(u[1] => 1 // 1, u[2] => 0 // 1),))

            @test isequal(terms[1],
                          1 // 48 * dt^5 + 31 // 640 * dt^7 + 8969 // 80640 * dt^9)
            @test isequal(terms[2],
                          1 + 1 // 12 * dt^2 + 1 // 20 * dt^4 + 127 // 2016 * dt^6 +
                          1513 // 12960 * dt^8)
        end
    end # @testset "nonlinear oscillator"

    @testset "integer coefficients" begin
        # reported by David Ketcheson on 2021-12-09
        A = reshape([1], 1, 1)
        b = [1]
        c = [1]
        @inferred modified_equation(A, b, c, 4)
    end

    @testset "average vector field method" begin
        series = bseries(3) do t, series
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

        diff = series - ExactSolution(series)
        @test mapreduce(abs ∘ last, +, diff) == 1 // 12
        @test @inferred(order_of_accuracy(series)) == 2
    end

    @testset "Runge-Kutta methods interface" begin
        # classical RK4
        A = [0 0 0 0
             1//2 0 0 0
             0 1//2 0 0
             0 0 1 0]
        b = [1 // 6, 1 // 3, 1 // 3, 1 // 6]
        rk = RungeKuttaMethod(A, b)

        # fourth-order accurate
        series_integrator = @inferred bseries(rk, 4)
        series_exact = @inferred ExactSolution(series_integrator)
        @test mapreduce(==, &, values(series_integrator), values(series_exact))
        @test @inferred(order_of_accuracy(series_integrator)) == 4

        # not fifth-order accurate
        series_integrator = @inferred bseries(rk, 5)
        series_exact = @inferred ExactSolution(series_integrator)
        @test mapreduce(==, &, values(series_integrator), values(series_exact)) == false
        @test @inferred(order_of_accuracy(series_integrator)) == 4
    end # @testset "Runge-Kutta methods interface"

    @testset "additive Runge-Kutta methods interface" begin
        @testset "IMEX Euler (partitioned Euler, symplectic Euler)" begin
            ex_euler = RungeKuttaMethod(@SMatrix([0]), @SVector [1])
            im_euler = RungeKuttaMethod(@SMatrix([1]), @SVector [1])
            ark = AdditiveRungeKuttaMethod([im_euler, ex_euler])

            # first-order accurate
            series_integrator = @inferred bseries(ark, 1)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(==, &, values(series_integrator), values(series_exact))
            @test @inferred(order_of_accuracy(series_integrator)) == 1

            # not second-order accurate
            series_integrator = @inferred bseries(ark, 2)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(==, &, values(series_integrator), values(series_exact)) == false

            # modified equations and modifying integrators
            let series_integrator = @inferred bseries(ark, 3)
                @testset "modified_equation" begin
                    mod_eq = @inferred modified_equation(series_integrator)

                    # Hairer, Lubich, Wanner (2006) Geometric numerical integration
                    # Table IX.10.1, p. 383
                    # Black nodes are `0`, white nodes are `1`
                    mod_eq_reference = Dict(rootedtree(Int[], Bool[]) => 0 // 1,
                                            rootedtree([1], Bool[0]) => 1 // 1,
                                            rootedtree([1], Bool[1]) => 1 // 1,
                                            rootedtree([1, 2], Bool[0, 0]) => 1 // 2,
                                            rootedtree([1, 2], Bool[1, 0]) => 1 // 2,
                                            rootedtree([1, 2], Bool[0, 1]) => -1 // 2,
                                            rootedtree([1, 2], Bool[1, 1]) => -1 // 2,
                                            rootedtree([1, 2, 3], Bool[0, 0, 0]) => 1 // 3,
                                            rootedtree([1, 2, 3], Bool[1, 0, 0]) => 1 // 3,
                                            rootedtree([1, 2, 3], Bool[0, 1, 0]) => -1 // 6,
                                            rootedtree([1, 2, 3], Bool[1, 1, 0]) => -1 // 6,
                                            rootedtree([1, 2, 3], Bool[0, 0, 1]) => -1 // 6,
                                            rootedtree([1, 2, 3], Bool[1, 0, 1]) => -1 // 6,
                                            rootedtree([1, 2, 3], Bool[0, 1, 1]) => 1 // 3,
                                            rootedtree([1, 2, 3], Bool[1, 1, 1]) => 1 // 3,
                                            rootedtree([1, 2, 2], Bool[0, 0, 0]) => 1 // 6,
                                            rootedtree([1, 2, 2], Bool[1, 0, 0]) => 1 // 6,
                                            rootedtree([1, 2, 2], Bool[0, 1, 0]) => -1 // 3,
                                            rootedtree([1, 2, 2], Bool[1, 1, 0]) => -1 // 3,
                                            rootedtree([1, 2, 2], Bool[0, 1, 1]) => 1 // 6,
                                            rootedtree([1, 2, 2], Bool[1, 1, 1]) => 1 // 6)
                    for (key, val) in mod_eq
                        @test val == mod_eq_reference[key]
                    end
                end

                @testset "modifying_integrator" begin
                    mod_int = @inferred modifying_integrator(series_integrator)
                end
            end
        end

        @testset "Störmer-Verlet" begin
            # Hairer, Lubich, Wanner (2002)
            # Geometric numerical integration
            # Table II.2.1
            As = [
                [0 0; 1//2 1//2],
                [1//2 0; 1//2 0],
            ]
            bs = [
                [1 // 2, 1 // 2],
                [1 // 2, 1 // 2],
            ]
            ark = AdditiveRungeKuttaMethod(As, bs)

            # second-order accurate
            series_integrator = @inferred bseries(ark, 2)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(==, &, values(series_integrator), values(series_exact))
            @test @inferred(order_of_accuracy(series_integrator)) == 2

            # not third-order accurate
            series_integrator = @inferred bseries(ark, 3)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(==, &, values(series_integrator), values(series_exact)) == false

            # modified equations and modifying integrators
            let series_integrator = @inferred bseries(ark, 4)
                @testset "modified_equation" begin
                    mod_eq = @inferred modified_equation(series_integrator)

                    # second-order accurate
                    for t in BicoloredRootedTreeIterator(2)
                        @test iszero(mod_eq[t])
                    end

                    # Hairer, Lubich, Wanner (2006) Geometric numerical integration
                    # Table IX.10.2, p. 386
                    # Black nodes are `0`, white nodes are `1`
                    series_reference = Dict(rootedtree(Int[], Bool[]) => 1 // 1,
                                            rootedtree([1], Bool[0]) => 1 // 1,
                                            rootedtree([1], Bool[1]) => 1 // 1,
                                            rootedtree([1, 2], Bool[0, 0]) => 1 // 2,
                                            rootedtree([1, 2], Bool[1, 0]) => 1 // 2,
                                            rootedtree([1, 2], Bool[0, 1]) => 1 // 2,
                                            rootedtree([1, 2], Bool[1, 1]) => 1 // 2,
                                            rootedtree([1, 2, 3], Bool[0, 0, 0]) => 1 // 4,
                                            rootedtree([1, 2, 3], Bool[1, 0, 0]) => 1 // 4,
                                            rootedtree([1, 2, 3], Bool[0, 1, 0]) => 0 // 1,
                                            rootedtree([1, 2, 3], Bool[1, 1, 0]) => 0 // 1,
                                            rootedtree([1, 2, 3], Bool[0, 0, 1]) => 1 // 4,
                                            rootedtree([1, 2, 3], Bool[1, 0, 1]) => 1 // 4,
                                            rootedtree([1, 2, 3], Bool[0, 1, 1]) => 1 // 4,
                                            rootedtree([1, 2, 3], Bool[1, 1, 1]) => 1 // 4,
                                            rootedtree([1, 2, 2], Bool[0, 0, 0]) => 1 // 2,
                                            rootedtree([1, 2, 2], Bool[1, 0, 0]) => 1 // 2,
                                            rootedtree([1, 2, 2], Bool[0, 1, 0]) => 1 // 4,
                                            rootedtree([1, 2, 2], Bool[1, 1, 0]) => 1 // 4,
                                            rootedtree([1, 2, 2], Bool[0, 1, 1]) => 1 // 4,
                                            rootedtree([1, 2, 2], Bool[1, 1, 1]) => 1 // 4)
                    for t in keys(series_reference)
                        @test series_reference[t] == series_integrator[t]
                    end
                    mod_eq_reference = Dict(rootedtree(Int[], Bool[]) => 0 // 1,
                                            rootedtree([1], Bool[0]) => 1 // 1,
                                            rootedtree([1], Bool[1]) => 1 // 1,
                                            rootedtree([1, 2], Bool[0, 0]) => 0 // 1,
                                            rootedtree([1, 2], Bool[1, 0]) => 0 // 1,
                                            rootedtree([1, 2], Bool[0, 1]) => 0 // 1,
                                            rootedtree([1, 2], Bool[1, 1]) => 0 // 1,
                                            rootedtree([1, 2, 3], Bool[0, 0, 0]) => 1 // 12,
                                            rootedtree([1, 2, 3], Bool[1, 0, 0]) => 1 // 12,
                                            rootedtree([1, 2, 3], Bool[0, 1, 0]) => -1 // 6,
                                            rootedtree([1, 2, 3], Bool[1, 1, 0]) => -1 // 6,
                                            rootedtree([1, 2, 3], Bool[0, 0, 1]) => 1 // 12,
                                            rootedtree([1, 2, 3], Bool[1, 0, 1]) => 1 // 12,
                                            rootedtree([1, 2, 3], Bool[0, 1, 1]) => 1 // 12,
                                            rootedtree([1, 2, 3], Bool[1, 1, 1]) => 1 // 12,
                                            rootedtree([1, 2, 2], Bool[0, 0, 0]) => 1 // 6,
                                            rootedtree([1, 2, 2], Bool[1, 0, 0]) => 1 // 6,
                                            rootedtree([1, 2, 2], Bool[0, 1, 0]) => -1 // 12,
                                            rootedtree([1, 2, 2], Bool[1, 1, 0]) => -1 // 12,
                                            rootedtree([1, 2, 2], Bool[0, 1, 1]) => -1 // 12,
                                            rootedtree([1, 2, 2], Bool[1, 1, 1]) => -1 // 12)
                    for t in keys(mod_eq_reference)
                        @test mod_eq_reference[t] == mod_eq[t]
                    end

                    # Hairer, Lubich, Wanner (2003)
                    # Geometric numerical integration illustrated by the Störmer-Verlet method
                    # https://doi.org/10.1017/S0962492902000144
                    # equation (4.8)
                    @testset "Pendulum, SymEngine" begin
                        dt = SymEngine.symbols("dt")
                        q, v = u = SymEngine.symbols("q, v")
                        f = ([v, 0], [0, -sin(q)])

                        series = modified_equation(f, u, dt, series_integrator)
                        series_reference = f[1] + f[2] +
                                           1 // 12 * dt^2 *
                                           [2 * cos(q) * v, sin(q) * cos(q) + sin(q) * v^2]
                        @test mapreduce(iszero ∘ SymEngine.expand, &,
                                        series - series_reference)
                    end
                end # @testset "modified_equation"

                @testset "modifying_integrator" begin
                    mod_int = @inferred modifying_integrator(series_integrator)

                    for t in BicoloredRootedTreeIterator(2)
                        @test iszero(mod_int[t])
                    end
                end
            end
        end

        @testset "SSPRK33 two times" begin
            # Using the same method multiple times is equivalent to using a plain RK
            # method without any splitting/partitioning/decomposition
            A = @SArray [0 0 0; 1 0 0; 1//4 1//4 0]
            b = @SArray [1 // 6, 1 // 6, 2 // 3]
            rk = RungeKuttaMethod(A, b)
            ark = AdditiveRungeKuttaMethod([rk, rk])

            # third-order accurate
            series_integrator = @inferred bseries(ark, 3)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(==, &, values(series_integrator), values(series_exact))
            @test @inferred(order_of_accuracy(series_integrator)) == 3

            # not fourth-order accurate
            series_integrator = @inferred bseries(ark, 4)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(==, &, values(series_integrator), values(series_exact)) == false

            # modified equations and modifying integrators
            let series_rk = bseries(rk, 5)
                series_ark = bseries(ark, order(series_rk))

                for colored_tree in keys(series_ark)
                    tree = rootedtree(colored_tree.level_sequence)
                    @test series_rk[tree] == series_ark[colored_tree]
                end

                @testset "modified_equation" begin
                    mod_eq_rk = modified_equation(series_rk)
                    mod_eq_ark = modified_equation(series_ark)
                    for colored_tree in keys(mod_eq_ark)
                        tree = rootedtree(colored_tree.level_sequence)
                        @test mod_eq_rk[tree] == mod_eq_ark[colored_tree]
                    end
                end

                @testset "modifying_integrator" begin
                    mod_int_rk = modifying_integrator(series_rk)
                    mod_int_ark = modifying_integrator(series_ark)
                    for colored_tree in keys(mod_int_ark)
                        tree = rootedtree(colored_tree.level_sequence)
                        @test mod_int_rk[tree] == mod_int_ark[colored_tree]
                    end
                end
            end
        end
    end # @testset "additive Runge-Kutta methods interface"

    @testset "Rosenbrock methods interface" begin
        @testset "Kaps, Rentrop (1979): GRK4A" begin
            # Kaps, Rentrop (1979)
            # Generalized Runge-Kutta methods of order four with stepsize control
            # for stiff ordinary differential equations
            # https://doi.org/10.1007/BF01396495
            Γ = [0.395 0 0 0;
                 -0.767672395484 0.395 0 0;
                 -0.851675323742 0.522967289188 0.395 0;
                 0.288463109545 0.880214273381e-1 -0.337389840627 0.395]
            A = [0 0 0 0;
                 0.438 0 0 0;
                 0.796920457938 0.730795420615e-1 0 0;
                 0.796920457938 0.730795420615e-1 0 0]
            b = [0.199293275701, 0.482645235674, 0.680614886256e-1, 0.25]
            ros = @inferred RosenbrockMethod(Γ, A, b)

            # fourth-order accurate
            series_integrator = @inferred bseries(ros, 5)
            @test @inferred(order_of_accuracy(series_integrator)) == 4

            # not fifth-order accurate
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(isapprox, &, values(series_integrator), values(series_exact)) ==
                  false
        end

        @testset "van Veldhuizen (1984)" begin
            # van Veldhuizen (1984)
            # D-stability and Kaps-Rentrop methods
            # https://doi.org/10.1007/BF02243574
            # Γ = [1//2 0 0 0;
            #      -4 1//2 0 0;
            #      -4 -1//2 1//2 0;
            #      1//4 -1//4 1 1//2]
            # A = [0 0 0 0;
            #      1 0 0 0;
            #      7//8 1//8 0 0;
            #      7//8 1//8 0 0]
            # b = [4//6, 2//6, -4//6, 4//6]
            # ros = @inferred RosenbrockMethod(Γ, A, b)
            # However, this does not work directly. Thus, we reverse-engineer
            # the coefficients as follows.
            #
            # Hairer, Wanner
            # Solving ODEs II
            # Implementation of Rosenbrock-type methods in Section IV.7.
            # The coefficients are transformed as
            # - C = I / γ - inv(Γ)
            # - A = A / Γ
            # - b' = b' / Γ
            # to yield the coefficients used in
            # http://www.unige.ch/~hairer/prog/stiff/Oldies/ros4.f
            C = [0 0 0 0;
                 -8 0 0 0;
                 -8 -1 0 0;
                 1//2 -1//2 2 0]
            γ = 1 // 2
            Γ = inv(I / γ - C)
            A = [0 0 0 0;
                 2 0 0 0;
                 7//4 1//4 0 0;
                 7//4 1//4 0 0]
            A = A * Γ
            b = [4 // 3, 2 // 3, -4 // 3, 4 // 3]
            b = (b' * Γ)'
            ros = @inferred RosenbrockMethod(Γ, A, b)

            # fourth-order accurate
            series_integrator = @inferred bseries(ros, 5)
            @test @inferred(order_of_accuracy(series_integrator)) == 4

            # not fifth-order accurate
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(isapprox, &, values(series_integrator), values(series_exact)) ==
                  false
        end
    end # @testset "Rosenbrock methods interface"

    @testset "Continuous stage Runge-Kutta methods interface" begin
        @testset "Average vector field method" begin
            M = fill(1 // 1, 1, 1)
            csrk = @inferred ContinuousStageRungeKuttaMethod(M)

            # TODO: This is no type stable at the moment
            # series = @inferred bseries(csrk, 8)
            series = bseries(csrk, 8)
            series_avf = @inferred bseries(AverageVectorFieldMethod(), order(series))
            @test all(iszero, values(series - series_avf))
        end

        @testset "Average vector field method with integer coefficient" begin
            M = fill(1, 1, 1)
            csrk = @inferred ContinuousStageRungeKuttaMethod(M)

            # TODO: This is no type stable at the moment
            # series = @inferred bseries(csrk, 8)
            series = bseries(csrk, 8)
            series_avf = @inferred bseries(AverageVectorFieldMethod(), order(series))
            @test all(iszero, values(series - series_avf))
        end

        @testset "Example in Section 4.2 of Miyatake and Butcher (2016)" begin
            # - Yuto Miyatake and John C. Butcher.
            #   "A characterization of energy-preserving methods and the construction of
            #   parallel integrators for Hamiltonian systems."
            #   SIAM Journal on Numerical Analysis 54, no. 3 (2016):
            #   [DOI: 10.1137/15M1020861](https://doi.org/10.1137/15M1020861)
            M = [-6//5 72//5 -36//1 24//1
                 72//5 -144//5 -48//1 72//1
                 -36//1 -48//1 720//1 -720//1
                 24//1 72//1 -720//1 720//1]
            csrk = @inferred ContinuousStageRungeKuttaMethod(M)

            # TODO: This is no type stable at the moment
            # series = @inferred bseries(csrk, 6)
            series = bseries(csrk, 6)

            @test @inferred(order_of_accuracy(series)) == 4
            @test is_energy_preserving(series)

            # Now with floating point coefficients
            M = Float64.(M)
            csrk = @inferred ContinuousStageRungeKuttaMethod(M)

            # TODO: This is no type stable at the moment
            # series = @inferred bseries(csrk, 6)
            series = bseries(csrk, 6)

            @test @inferred(order_of_accuracy(series)) == 4
            @test_broken is_energy_preserving(series)
        end

        @testset "SymEngine.jl" begin
            # Examples in Section 5.3.1
            α = SymEngine.symbols("α")
            α1 = 1 / (36 * α - 7)
            M = [α1+4 -6 * α1-6 6*α1
                 -6 * α1-6 36 * α1+12 -36*α1
                 6*α1 -36*α1 36*α1]
            csrk = @inferred ContinuousStageRungeKuttaMethod(M)

            # TODO: This is no type stable at the moment
            # series = @inferred bseries(csrk, 5)
            series = bseries(csrk, 5)

            # The simple test does not work at the moment due to missing
            # simplifications in SymEngine.jl
            @test_broken @inferred(order_of_accuracy(series)) == 4
            exact = @inferred ExactSolution(series)
            for o in 1:4
                for t in RootedTreeIterator(o)
                    expr = SymEngine.expand(series[t] - exact[t])
                    @test iszero(SymEngine.expand(1 * expr))
                end
            end

            # TODO: This is currently not implemented
            @test_broken is_energy_preserving(series)
        end

        @testset "SymPy.jl" begin
            # Examples in Section 5.3.1
            α = SymPy.symbols("α", real = true)
            α1 = 1 / (36 * α - 7)
            M = [α1+4 -6 * α1-6 6*α1
                 -6 * α1-6 36 * α1+12 -36*α1
                 6*α1 -36*α1 36*α1]
            csrk = @inferred ContinuousStageRungeKuttaMethod(M)

            # TODO: This is no type stable at the moment
            # series = @inferred bseries(csrk, 5)
            series = bseries(csrk, 5)

            @test @inferred(order_of_accuracy(series)) == 4

            # TODO: This is currently not implemented
            @test_broken is_energy_preserving(series)
        end

        @testset "Symbolics.jl" begin
            # Examples in Section 5.3.1
            Symbolics.@variables α
            α1 = 1 / (36 * α - 7)
            M = [α1+4 -6 * α1-6 6*α1
                 -6 * α1-6 36 * α1+12 -36*α1
                 6*α1 -36*α1 36*α1]
            csrk = @inferred ContinuousStageRungeKuttaMethod(M)

            # TODO: This is no type stable at the moment
            # series = @inferred bseries(csrk, 2)
            series = bseries(csrk, 2)

            # TODO: There are errors from Symbolics.jl if we go to higher
            #       orders.
            # @test @inferred(order_of_accuracy(series)) == 4

            # # TODO: This is currently not implemented
            @test_broken is_energy_preserving(series)
        end
    end

    @testset "multirate infinitesimal split methods interface" begin
        # NOTE: This is still considered experimental and might change at any time!

        @testset "KW43" begin
            # Oswald Knoth, Ralf Wolke
            # Implicit-explicit Runge-Kutta methods for computing atmospheric reactive flows
            # Applied Numerical Mathematics, Volume 28, (1998) Issue 2-4
            # https://doi.org/10.1016/S0168-9274(98)00051-8
            ns = 4
            A = zeros(Polynomial{Rational{Int}, :x}, ns + 1, ns)
            D = zeros(Int, ns + 1, ns)
            G = zeros(Rational{Int}, ns + 1, ns)
            c = zeros(Rational{Int}, ns + 1)
            dts = zeros(Rational{Int}, ns + 1)

            A[2, 1] = 1 // 2
            A[3, 1] = -1 // 6
            A[3, 2] = 2 // 3
            A[4, 1] = 1 // 3
            A[4, 2] = -1 // 3
            A[4, 3] = 1
            A[5, 1] = 1 // 6
            A[5, 2] = 1 // 3
            A[5, 3] = 1 // 3
            A[5, 4] = 1 // 6
            A[5, :] = A[5, :] - A[4, :]
            A[4, :] = A[4, :] - A[3, :]
            A[3, :] = A[3, :] - A[2, :]

            c[2] = 1 // 2
            c[3] = 1 // 2
            c[4] = 1
            c[5] = 1

            D[2, 1] = 1
            D[3, 2] = 1
            D[4, 3] = 1
            D[5, 4] = 1

            dts[1] = 0
            dts[2] = 1 // 2
            dts[3] = 0
            dts[4] = 1 // 2
            dts[5] = 0

            mis = @inferred MultirateInfinitesimalSplitMethod(A, D, G, c)

            # third-order accurate
            series_integrator = @inferred bseries(mis, 3)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(==, &, values(series_integrator), values(series_exact))

            # not fourth-order accurate
            series_integrator = @inferred bseries(mis, 4)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(==, &, values(series_integrator), values(series_exact)) == false
        end

        @testset "MRI-GARK-ERK33" begin
            # Sandu, Adrian.
            # "A class of multirate infinitesimal GARK methods".
            # SIAM Journal on Numerical Analysis 57, no. 5 (2019): 2300-2327.
            # https://doi.org/10.1137/18M1205492
            ns = 3
            A = zeros(Polynomial{Float64, :x}, ns + 1, ns)
            D = zeros(ns + 1, ns)
            G = zeros(ns + 1, ns)
            c = zeros(ns + 1)
            dts = zeros(ns + 1)

            delta = -0.5

            # Γ^i in Sandu (2019) is the i-th term
            A[2, 1] = Polynomial([1 / 3])
            A[3, 1] = Polynomial([(-6 * delta - 7) / 12, (2 * delta + 1) / 2])
            A[3, 2] = Polynomial([(6 * delta + 11) / 12, -(2 * delta + 1) / 2])
            A[4, 1] = Polynomial([0.0, 0.5])
            A[4, 2] = Polynomial([(6 * delta - 5) / 12, -(2 * delta + 1) / 2])
            A[4, 3] = Polynomial([(3 - 2 * delta) / 4, delta])

            # These should be the usual `c` coefficients → c[3] = 2 / 3?
            c[2] = 1 / 3
            c[3] = 2 / 3
            c[4] = 1

            # All one???
            D[2, 1] = 1
            D[3, 2] = 1
            D[4, 3] = 1

            # Differences of `c`s?
            dts[1] = 0
            dts[2] = 1 / 3
            dts[3] = 1 / 3
            dts[4] = 1 / 3

            mis = @inferred MultirateInfinitesimalSplitMethod(A, D, G, c)

            # third-order accurate
            series_integrator = @inferred bseries(mis, 3)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(≈, &, values(series_integrator), values(series_exact))

            # not fourth-order accurate
            series_integrator = @inferred bseries(mis, 4)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(≈, &, values(series_integrator), values(series_exact)) == false
        end

        @testset "MRI-GARK-ERK45a" begin
            # Sandu, Adrian.
            # "A class of multirate infinitesimal GARK methods".
            # SIAM Journal on Numerical Analysis 57, no. 5 (2019): 2300-2327.
            # https://doi.org/10.1137/18M1205492
            ns = 5
            A = zeros(Polynomial{Float64, :x}, ns + 1, ns)
            D = zeros(ns + 1, ns)
            G = zeros(ns + 1, ns)
            c = zeros(ns + 1)
            dts = zeros(ns + 1)

            # Γ^i in Sandu (2019) is the i-th term
            A[2, 1] = Polynomial([1 / 5])

            A[3, 1] = Polynomial([-53 / 16, 503 / 80])
            A[3, 2] = Polynomial([281 / 80, -503 / 80])

            A[4, 1] = Polynomial([-36562993.0 / 71394880.0, -1365537.0 / 35697440.0])
            A[4, 2] = Polynomial([34903117.0 / 17848720.0, 4963773.0 / 7139488.0])
            A[4, 3] = Polynomial([-88770499.0 / 71394880.0, -1465833.0 / 2231090.0])

            A[5, 1] = Polynomial([-7631593.0 / 71394880.0, 66974357.0 / 35697440.0])
            A[5, 2] = Polynomial([-166232021.0 / 35697440.0, 21445367.0 / 7139488.0])
            A[5, 3] = Polynomial([6068517.0 / 1519040.0, -3.0])
            A[5, 4] = Polynomial([8644289.0 / 8924360.0, -8388609.0 / 4462180.0])

            A[6, 1] = Polynomial([277061.0 / 303808.0, -18227.0 / 7520.0])
            A[6, 2] = Polynomial([-209323.0 / 1139280.0, 2.0])
            A[6, 3] = Polynomial([-1360217.0 / 1139280.0, 1.0])
            A[6, 4] = Polynomial([-148789.0 / 56964.0, 5.0])
            A[6, 5] = Polynomial([147889.0 / 45120.0, -41933.0 / 7520.0])

            # These should be the usual `c` coefficients
            c[2] = 1 / 5
            c[3] = 2 / 5
            c[4] = 3 / 5
            c[5] = 4 / 5
            c[6] = 1

            D[2, 1] = 0
            D[3, 2] = 1
            D[4, 3] = 1
            D[5, 4] = 1
            D[6, 5] = 1

            # Differences of `c`s?
            dts[1] = 0.0
            dts[2] = 1 / 5
            dts[3] = 1 / 5
            dts[4] = 1 / 5
            dts[5] = 1 / 5
            dts[6] = 1 / 5

            mis = @inferred MultirateInfinitesimalSplitMethod(A, D, G, c)

            # fourth-order accurate
            series_integrator = @inferred bseries(mis, 4)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(≈, &, values(series_integrator), values(series_exact))

            # not fifth-order accurate
            series_integrator = @inferred bseries(mis, 5)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(≈, &, values(series_integrator), values(series_exact)) == false
        end

        @testset "EB4" begin
            # Krogstad, Stein
            # "Generalized integrating factor methods for stiff PDEs".
            # Journal of Computational Physics 203, no. 1 (2005): 72-88.
            # https://doi.org/10.1016/j.jcp.2004.08.006
            ns = 4
            A = zeros(Polynomial{Float64, :x}, ns + 1, ns)
            D = zeros(ns + 1, ns)
            G = zeros(ns + 1, ns)
            c = zeros(ns + 1)
            dts = zeros(ns + 1)

            fac1 = 1.0
            fac2 = 2.0

            A[2, 1] = Polynomial([0.5])

            A[3, 1] = Polynomial([0.5, -1 / fac1])
            A[3, 2] = Polynomial([0.0, 1 / fac1])

            A[4, 1] = Polynomial([1.0, -2 / fac1])
            A[4, 3] = Polynomial([0.0, 2 / fac1])

            A[5, 1] = Polynomial([1.0, -3 / fac1, 4 / fac2])
            A[5, 2] = Polynomial([0.0, 2 / fac1, -4 / fac2])
            A[5, 3] = Polynomial([0.0, 2 / fac1, -4 / fac2])
            A[5, 4] = Polynomial([0.0, -1 / fac1, 4 / fac2])

            c[1] = 0
            c[2] = 0.5
            c[3] = 0.5
            c[4] = 1
            c[5] = 1

            mis = @inferred MultirateInfinitesimalSplitMethod(A, D, G, c)

            # fourth-order accurate
            series_integrator = @inferred bseries(mis, 4)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(≈, &, values(series_integrator), values(series_exact))

            # not fifth-order accurate
            series_integrator = @inferred bseries(mis, 5)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(≈, &, values(series_integrator), values(series_exact)) == false
        end

        @testset "OwrenRK4" begin
            # Owren, Brynjulf.
            # "Order conditions for commutator-free Lie group methods".
            # Journal of Physics A: Mathematical and General 39, no. 19 (2006): 5585.
            # https://doi.org/10.1088/0305-4470/39/19/S15
            ns = 5
            A = zeros(Polynomial{Float64, :x}, ns + 1, ns)
            D = zeros(ns + 1, ns)
            G = zeros(ns + 1, ns)
            c = zeros(ns + 1)
            dts = zeros(ns + 1)

            A[2, 1] = Polynomial([1 / 2])

            A[3, 2] = Polynomial([1 / 2])

            A[4, 1] = Polynomial([-1 / 2])
            A[4, 3] = Polynomial([1.0])

            A[5, 1] = Polynomial([1 / 4])
            A[5, 2] = Polynomial([1 / 6])
            A[5, 3] = Polynomial([1 / 6])
            A[5, 4] = Polynomial([-1 / 12])

            A[6, 1] = Polynomial([-1 / 12])
            A[6, 2] = Polynomial([1 / 6])
            A[6, 3] = Polynomial([1 / 6])
            A[6, 4] = Polynomial([1 / 4])

            D[4, 2] = 1
            D[6, 5] = 1

            dts[1] = 0
            dts[2] = 1 / 2
            dts[3] = 1 / 2
            dts[4] = 1 / 2
            dts[5] = 1 / 2
            dts[6] = 1 / 2

            mis = @inferred MultirateInfinitesimalSplitMethod(A, D, G, c)

            # fourth-order accurate
            series_integrator = @inferred bseries(mis, 4)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(≈, &, values(series_integrator), values(series_exact))

            # not fifth-order accurate
            series_integrator = @inferred bseries(mis, 5)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(≈, &, values(series_integrator), values(series_exact)) == false
        end

        @testset "MISBI4_4" begin
            # Knoth, Oswald, and Joerg Wensch.
            # "Generalized split-explicit Runge–Kutta methods for the compressible Euler equations".
            # Monthly Weather Review 142, no. 5 (2014): 2067-2081.
            # https://doi.org/10.1175/MWR-D-13-00068.1
            # Table 4, method MIS4a
            ns = 4
            A = zeros(Polynomial{Float64, :x}, ns + 1, ns)
            D = zeros(ns + 1, ns)
            G = zeros(ns + 1, ns)
            c = zeros(ns + 1)
            dts = zeros(ns + 1)

            # β in Knoth, Wensch (2014)
            A[2, 1] = Polynomial([1 / 2])
            A[2, 1] = Polynomial([0.38758444641450318])
            A[3, 1] = Polynomial([-2.5318448354142823e-002])
            A[3, 2] = Polynomial([0.38668943087310403])
            A[4, 1] = Polynomial([0.20899983523553325])
            A[4, 2] = Polynomial([-0.45856648476371231])
            A[4, 3] = Polynomial([0.43423187573425748])
            A[5, 1] = Polynomial([-0.10048822195663100])
            A[5, 2] = Polynomial([-0.46186171956333327])
            A[5, 3] = Polynomial([0.83045062122462809])
            A[5, 4] = Polynomial([0.27014914900250392])

            # seems to be unused?
            c[2] = 0.38758444641450318
            c[3] = 0.61521685655017821
            c[4] = 0.23254717315441453
            c[5] = 1.0000000000000002

            # α in Knoth, Wensch (2014)
            D[3, 2] = 0.52349249922385610
            D[4, 2] = 1.1683374366893629
            D[4, 3] = -0.75762080241712637
            D[5, 2] = -3.6477233846797109e-002
            D[5, 3] = 0.56936148730740477
            D[5, 4] = 0.47746263002599681

            # γ in Knoth, Wensch (2014)
            G[3, 2] = 0.13145089796226542
            G[4, 2] = -0.36855857648747881
            G[4, 3] = 0.33159232636600550
            G[5, 2] = -6.5767130537473045e-002
            G[5, 3] = 4.0591093109036858e-002
            G[5, 4] = 6.4902111640806712e-002

            for i in 1:(ns + 1)
                for j in 1:(i - 1)
                    dts[i] = dts[i] + A[i, j][0]
                end
            end

            mis = @inferred MultirateInfinitesimalSplitMethod(A, D, G, c)

            # third-order accurate
            series_integrator = @inferred bseries(mis, 3)
            series_exact = @inferred ExactSolution(series_integrator)
            @test_broken mapreduce(≈, &, values(series_integrator), values(series_exact))

            # not fourth-order accurate
            series_integrator = @inferred bseries(mis, 4)
            series_exact = @inferred ExactSolution(series_integrator)
            @test mapreduce(≈, &, values(series_integrator), values(series_exact)) == false
        end
    end # @testset "multirate infinitesimal split methods interface"

    @testset "Energy preservation (Hamiltonian systems)" begin
        @testset "Pseudo-energy-preserving order 4" begin
            # References
            # Celledoni, Elena; McLachlan, Robert I.; McLaren, David I.;
            # Owren, Brynjulf; G. Reinout W. Quispel; Wright, William M.
            # Energy-preserving Runge-Kutta methods.
            # ESAIM: Mathematical Modelling and Numerical Analysis -
            # Modélisation Mathématique et Analyse Numérique,
            # Volume 43 (2009) no. 4, pp. 645-649.
            # doi : 10.1051/m2an/2009020. http://www.numdam.org/articles/10.1051/m2an/2009020/
            A = [0 0 0
                 1//3 0 0
                 -5//48 15//16 0]
            b = [1 // 10, 1 // 2, 2 // 5]
            rk = RungeKuttaMethod(A, b)

            # This method is E-P up to order 4
            @test energy_preserving_order(rk, 10) == 4
        end

        @testset "Pseudo-energy-preserving order 3" begin
            # References
            # Celledoni, Elena; McLachlan, Robert I.; McLaren, David I.;
            # Owren, Brynjulf; G. Reinout W. Quispel; Wright, William M.
            # Energy-preserving Runge-Kutta methods.
            # ESAIM: Mathematical Modelling and Numerical Analysis -
            # Modélisation Mathématique et Analyse Numérique,
            # Volume 43 (2009) no. 4, pp. 645-649.
            # doi : 10.1051/m2an/2009020. http://www.numdam.org/articles/10.1051/m2an/2009020/
            A = [0 0
                 2//3 0]
            b = [1 // 4, 3 // 4]
            rk = RungeKuttaMethod(A, b)

            # This method is E-P up to order 3
            @test energy_preserving_order(rk, 10) == 3
        end

        @testset "Classical RK Method" begin
            A = [0//1 0//1 0//1 0//1
                 1//2 0//1 0//1 0//1
                 0//1 1//2 0//1 0//1
                 0//1 0//1 1//1 0//1]
            b = [1 // 6, 1 // 3, 1 // 3, 1 // 6]
            rk = RungeKuttaMethod(A, b)
            # This method is E-P up to order 4
            @test energy_preserving_order(rk, 10) == 4
        end

        @testset "Effective Order" begin
            # References
            # Butcher, J.C. (1969). The effective order of Runge-Kutta methods.
            # In: Morris, J.L. (eds) Conference on the Numerical Solution of
            # Differential Equations. Lecture Notes in Mathematics, vol 109.
            # Springer, Berlin, Heidelberg. https://doi.org/10.1007/BFb0060019
            A = [0 0 0 0 0
                 1//5 0 0 0 0
                 0 2//5 0 0 0
                 3//16 0 5//16 0 0
                 1//4 0 -5//4 2 0]

            b = [1 // 6, 0, 0, 2 // 3, 1 // 6]
            rk = RungeKuttaMethod(A, b)
            #This method is E-P up to order 4
            @test energy_preserving_order(rk, 10) == 4
        end

        @testset "Average Vector Field (AVF)" begin
            # select order
            p = 7
            series = bseries(p) do t, series
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
            @test is_energy_preserving(series) == true
        end

        @testset "Floating point coefficients" begin
            # AVF method again with various types of coefficients
            series = bseries(AverageVectorFieldMethod(Float32), 7)
            @test is_energy_preserving(series)

            series = bseries(AverageVectorFieldMethod(Float64), 7)
            @test is_energy_preserving(series)

            series = bseries(AverageVectorFieldMethod(BigFloat), 7)
            # TODO: This test is currently broken and throws an error
            @test_broken is_energy_preserving(series)
        end

        @testset "Symbolic coefficients" begin
            @testset "SymEngine.jl" begin
                # This method is second-order accurate. Thus, it is
                # energy-preserving up to order two.
                α = SymEngine.symbols("α")
                A = [0 0; 1/(2 * α) 0]
                b = [1 - α, α]
                c = [0, 1 / (2 * α)]
                series_integrator = @inferred(bseries(A, b, c, 2))
                @test @inferred(order_of_accuracy(series_integrator)) == 2
                # TODO: This test is currently broken and throws an error
                @test_broken is_energy_preserving(series_integrator)
            end

            @testset "SymPy.jl" begin
                # This method is second-order accurate. Thus, it is
                # energy-preserving up to order two.
                α = SymPy.symbols("α", real = true)
                A = [0 0; 1/(2 * α) 0]
                b = [1 - α, α]
                c = [0, 1 / (2 * α)]
                series_integrator = @inferred(bseries(A, b, c, 2))
                @test @inferred(order_of_accuracy(series_integrator)) == 2
                # TODO: This test is currently broken and throws an error
                @test_broken is_energy_preserving(series_integrator)
            end

            @testset "Symbolics.jl" begin
                # This method is second-order accurate. Thus, it is
                # energy-preserving up to order two.
                Symbolics.@variables α
                A = [0 0; 1/(2 * α) 0]
                b = [1 - α, α]
                c = [0, 1 / (2 * α)]
                series_integrator = @inferred(bseries(A, b, c, 2))
                @test @inferred(order_of_accuracy(series_integrator)) == 2
                # TODO: This test is currently broken and throws an error
                @test_broken is_energy_preserving(series_integrator)
            end
        end
    end

    @testset "Aqua" begin
        if VERSION < v"1.7"
            Aqua.test_all(BSeries;
                          # We would like to check for ambiguities but cannot do so right now because
                          # of https://github.com/JuliaTesting/Aqua.jl/issues/79
                          # Thus, we do not test for ambiguities here but run an additional test
                          # below excluding ambiguity tests with Base.
                          ambiguities = false)
        else
            Aqua.test_all(BSeries;
                          ambiguities = (; exclude = [
                                             getindex, # https://github.com/stevengj/LaTeXStrings.jl/issues/61
                                         ]),
                          # Requires.jl is not loaded on new versions of Julia
                          stale_deps = (; ignore = [:Requires]))
        end

        # No Base and as extra test for the reason described above
        @testset "ambiguities" begin
            Aqua.test_ambiguities([BSeries])
        end
    end
end # @testset "BSeries"
