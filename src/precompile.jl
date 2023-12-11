# We cannot use `@precompile_all_calls` from SnoopPrecompile.jl
# directly on code since the internal buffers of RootedTrees.jl are not set up
# during precompilation, see also
# https://github.com/SciML/RootedTrees.jl/pull/92
#
# Thus, we use the older tools from SnoopCompile.jl to generate precompile
# statements. This is based on the following code:
# ```julia
# julia> using SnoopCompile, ProfileView; tinf = @snoopi_deep begin
#
#   using BSeries
#
#   A = [0 0 0 0; 1//2 0 0 0; 0 1//2 0 0; 0 0 1 0];
#   b = [1//6, 1//3, 1//3, 1//6];
#   c = [0, 1//2, 1//2, 1];
#   series = bseries(A, b, c, 5)
#   exact = ExactSolution(series)
#   series - exact
#   order_of_accuracy(series)
#
#   As = [
#       [0 0; 1//2 1//2],
#       [1//2 0; 1//2 0],
#   ]
#   bs = [
#       [1 // 2, 1 // 2],
#       [1 // 2, 1 // 2],
#   ]
#   ark = AdditiveRungeKuttaMethod(As, bs)
#   series = bseries(ark, 3)
#   series - ExactSolution(series)
#   order_of_accuracy(series)
#   modified_equation(series)
#   modifying_integrator(series)
#
#   γ = [0.395 0 0 0;
#        -0.767672395484 0.395 0 0;
#        -0.851675323742 0.522967289188 0.395 0;
#        0.288463109545 0.880214273381e-1 -0.337389840627 0.395]
#   A = [0 0 0 0;
#        0.438 0 0 0;
#        0.796920457938 0.730795420615e-1 0 0;
#        0.796920457938 0.730795420615e-1 0 0]
#   b = [0.199293275701, 0.482645235674, 0.680614886256e-1, 0.25]
#   ros = RosenbrockMethod(γ, A, b)
#   series = bseries(ros, 5)
#   series - ExactSolution(series)
#   modified_equation(series)
#   modifying_integrator(series)
#
#   series = bseries(5) do t, series
#       if order(t) in (0, 1)
#           return 1 // 1
#       else
#           v = 1 // 1
#           n = 0
#           for subtree in SubtreeIterator(t)
#               v *= series[subtree]
#               n += 1
#           end
#           return v / (n + 1)
#       end
#   end
#   series - ExactSolution(series)
#   order_of_accuracy(series)
#   modified_equation(series)
#   modifying_integrator(series)
#
#   A = [0 0 0 0 0;
#        1//5 0 0 0 0;
#        0 2//5 0 0 0;
#        3//16 0 5//16 0 0;
#        1//4 0 -5//4 2 0]
#   b = [1 // 6, 0, 0, 2 // 3, 1 // 6]
#   rk_a = RungeKuttaMethod(A, b)
#   series_a = bseries(rk_a, 6)
#   order_of_accuracy(series_a)
#   A = [0 0 0 0 0;
#        1//5 0 0 0 0;
#        0 2//5 0 0 0;
#        75//64 -9//4 117//64 0 0;
#        -37//36 7//3 -3//4 4//9 0]
#   b = [19 // 144, 0, 25 // 48, 2 // 9, 1 // 8]
#   rk_b = RungeKuttaMethod(A, b)
#   series_b = bseries(rk_b, 6)
#   order_of_accuracy(series_b)
#   A = [0 0 0 0 0;
#        1//5 0 0 0 0;
#        0 2//5 0 0 0;
#        161//192 -19//12 287//192 0 0;
#        -27//28 19//7 -291//196 36//49 0]
#   b = [7 // 48, 0, 475 // 1008, 2 // 7, 7 // 72]
#   rk_c = RungeKuttaMethod(A, b)
#   series_c = bseries(rk_c, 6)
#   order_of_accuracy(series_c)
#   series = compose(series_b, series_a, series_c, normalize_stepsize = true)
#   order_of_accuracy(series)
#   modified_equation(series)
#   modifying_integrator(series)
#
#   end
#
# InferenceTimingNode: 3.093312/4.309741 on Core.Compiler.Timings.ROOT() with 63 direct children
#
# julia> ttot, pcs = SnoopCompile.parcel(tinf);
#
# julia> ttot
# 1.2164289380000002
#
# julia> SnoopCompile.write("/tmp/precompiles_BSeries", pcs, has_bodyfunction = true)
# Core: no precompile statements out of 0.000186417
# Base.Threads: precompiled 0.005340074 out of 0.005340074
# RootedTrees: precompiled 0.12003955399999999 out of 0.12143691799999999
# Base: precompiled 0.17134280800000004 out of 0.17531468100000003
# BSeries: precompiled 0.9096680570000002 out of 0.9097387219999999
# ```
# See https://timholy.github.io/SnoopCompile.jl/dev/snoopi_deep_parcel/

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(bseries), Matrix{Rational{Int64}}, Vector{Rational{Int64}},
                          Vector{Rational{Int64}}, Int64})   # time: 0.28764844
    Base.precompile(Tuple{typeof(bseries),
                          AdditiveRungeKuttaMethod{Rational{Int64},
                                                   Vector{RungeKuttaMethod{Rational{Int64},
                                                                           Matrix{Rational{Int64}},
                                                                           Vector{Rational{Int64}}}}},
                          Int64})   # time: 0.24456108
    Base.precompile(Tuple{typeof(bseries),
                          RosenbrockMethod{Float64, Matrix{Float64}, Vector{Float64}},
                          Int64})   # time: 0.11124245
    Base.precompile(Tuple{Core.kwftype(typeof(compose)),
                          NamedTuple{(:normalize_stepsize,), Tuple{Bool}}, typeof(compose),
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}},
                                           Rational{Int64}},
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}},
                                           Rational{Int64}},
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}},
                                           Rational{Int64}}})   # time: 0.052562617
    Base.precompile(Tuple{typeof(modified_equation),
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Float64}})   # time: 0.029089022
    Base.precompile(Tuple{typeof(modified_equation),
                          TruncatedBSeries{BicoloredRootedTree{Int64, Vector{Int64},
                                                               Vector{Bool}},
                                           Rational{Int64}}})   # time: 0.026970008
    Base.precompile(Tuple{typeof(-),
                          TruncatedBSeries{BicoloredRootedTree{Int64, Vector{Int64},
                                                               Vector{Bool}},
                                           Rational{Int64}},
                          TruncatedBSeries{BicoloredRootedTree{Int64, Vector{Int64},
                                                               Vector{Bool}},
                                           Rational{Int64}}})   # time: 0.024954015
    Base.precompile(Tuple{typeof(bseries), Function, Int64})   # time: 0.021058364
    Base.precompile(Tuple{typeof(order_of_accuracy),
                          TruncatedBSeries{BicoloredRootedTree{Int64, Vector{Int64},
                                                               Vector{Bool}},
                                           Rational{Int64}}})   # time: 0.018548366
    Base.precompile(Tuple{typeof(-),
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}},
                                           Rational{Int64}},
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}},
                                           Rational{Int64}}})   # time: 0.015549354
    Base.precompile(Tuple{Type{ExactSolution},
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}},
                                           Rational{Int64}}})   # time: 0.014666914
    Base.precompile(Tuple{typeof(order_of_accuracy),
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}},
                                           Rational{Int64}}})   # time: 0.012643992
    Base.precompile(Tuple{typeof(-),
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Float64},
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Float64}})   # time: 0.010377383
    Base.precompile(Tuple{typeof(modified_equation),
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}},
                                           Rational{Int64}}})   # time: 0.009575064
    Base.precompile(Tuple{Type{ExactSolution},
                          TruncatedBSeries{BicoloredRootedTree{Int64, Vector{Int64},
                                                               Vector{Bool}},
                                           Rational{Int64}}})   # time: 0.008866523
    Base.precompile(Tuple{typeof(getindex),
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}},
                                           Rational{Int64}},
                          RootedTree{Int64,
                                     SubArray{Int64, 1, Vector{Int64},
                                              Tuple{UnitRange{Int64}}, true}}})   # time: 0.008717881
    Base.precompile(Tuple{Type{ExactSolution},
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Float64}})   # time: 0.007198428
    Base.precompile(Tuple{typeof(modifying_integrator),
                          TruncatedBSeries{BicoloredRootedTree{Int64, Vector{Int64},
                                                               Vector{Bool}},
                                           Rational{Int64}}})   # time: 0.001966412
    Base.precompile(Tuple{typeof(modifying_integrator),
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}},
                                           Rational{Int64}}})   # time: 0.001889854
    Base.precompile(Tuple{typeof(modifying_integrator),
                          TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Float64}})   # time: 0.001581901
end

_precompile_()
