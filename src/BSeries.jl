module BSeries

@doc read(joinpath(dirname(@__DIR__), "README.md"), String) BSeries

import Base: +, -, *, /, \

using Reexport: @reexport

@reexport using RootedTrees
using RootedTrees: AbstractRootedTree

@reexport using OrderedCollections: OrderedDict

if !isdefined(Base, :get_extension)
    using Requires: @require
end

using Latexify: Latexify, LaTeXString
using Combinatorics: Combinatorics, permutations
using LinearAlgebra: LinearAlgebra, rank
using SparseArrays: SparseArrays, sparse

@reexport using Polynomials: Polynomials, Polynomial

export TruncatedBSeries, ExactSolution

export order_of_accuracy

export bseries, substitute, compose, evaluate

export modified_equation, modifying_integrator

export elementary_differentials

export AverageVectorFieldMethod

export ContinuousStageRungeKuttaMethod

export MultirateInfinitesimalSplitMethod

export renormalize!

export is_energy_preserving, energy_preserving_order

# Types used for traits
# These traits may decide between different algorithms based on the
# corresponding complexity etc.
abstract type AbstractEvaluation end
struct LazyEvaluation <: AbstractEvaluation end # lazy evaluation up to arbitrary order
struct EagerEvaluation <: AbstractEvaluation end # eager evaluation up to a fixed order
struct MemoizedEvaluation <: AbstractEvaluation end # lazy evaluation memoizing the results

evaluation_type(::Any) = LazyEvaluation() # assume lazy evaluation by default
evaluation_type(::AbstractDict) = EagerEvaluation() # dictionaries store results eagerly

# internal interface
# TODO: `bseries(series::AbstractBSeries) = series` or a `copy`?
# TODO: `bseries(series, order)` returns a `copy`? up to `order`, if available
# TODO: `truncate!(series, order)`

"""
    TruncatedBSeries

A struct that can describe B-series of both numerical integration methods
(where the coefficient of the empty tree is unity) and right-hand sides of an
ordinary differential equation and perturbations thereof (where the coefficient
of the empty tree is zero) up to a prescribed [`order`](@ref).

Generally, this kind of `struct` should be constructed via [`bseries`](@ref)
or one of the other functions returning a B-series, e.g.,
[`modified_equation`](@ref) or [`modifying_integrator`](@ref).
"""
struct TruncatedBSeries{T <: AbstractRootedTree, V} <: AbstractDict{T, V}
    coef::OrderedDict{T, V}
end

TruncatedBSeries{T, V}() where {T, V} = TruncatedBSeries{T, V}(OrderedDict{T, V}())

# general interface methods of `AbstractDict` for `TruncatedBSeries`
@inline Base.iterate(series::TruncatedBSeries) = iterate(series.coef)
@inline Base.iterate(series::TruncatedBSeries, state) = iterate(series.coef, state)

@inline Base.empty(series::TruncatedBSeries) = TruncatedBSeries(empty(series.coef))
@inline function Base.empty(series::TruncatedBSeries, ::Type{T}, ::Type{V}) where {T, V}
    TruncatedBSeries(empty(series.coef, T, V))
end
@inline Base.empty!(series::TruncatedBSeries) = (empty!(series.coef); series)

@inline Base.length(series::TruncatedBSeries) = length(series.coef)

@inline function Base.getindex(series::TruncatedBSeries, t::AbstractRootedTree)
    getindex(series.coef, t)
end
@inline function Base.setindex!(series::TruncatedBSeries, val, t::AbstractRootedTree)
    setindex!(series.coef, val, t)
end

@inline function Base.get(series::TruncatedBSeries, t::AbstractRootedTree, default)
    get(series.coef, t, default)
end
@inline function Base.get(f::F, series::TruncatedBSeries,
                          t::AbstractRootedTree) where {F <: Function}
    get(f, series.coef, t)
end

@inline function Base.getkey(series::TruncatedBSeries, t::AbstractRootedTree, default)
    getkey(series.coef, t, default)
end

@inline function Base.delete!(series::TruncatedBSeries, t::AbstractRootedTree)
    delete!(series.coef, t)
    return series
end

@inline Base.pop!(series::TruncatedBSeries, t::AbstractRootedTree) = pop!(series.coef, t)
@inline function Base.pop!(series::TruncatedBSeries, t::AbstractRootedTree, default)
    pop!(series.coef, t, default)
end

@inline Base.sizehint!(series::TruncatedBSeries, n) = sizehint!(series.coef, n)

# internal interface of B-series
# evaluation_type(::TruncatedBSeries) = EagerEvaluation() # is already the default for dicts

# TODO: TruncatedBSeries
# This assumes that the `coef` are always stored with increasing `order` of the
# rooted trees and that the underlying data format in OrderedCollections.jl does
# not change. A general fallback would be
#   RootedTrees.order(series::TruncatedBSeries) = maximum(order, keys(series.coef))
# but that is O(n) instead of O(1). Since we do not consider the constructor as
# public API but as internal implementation detail, users violating this assumption
# are outside of the public API and may run into self-made problems.
"""
    order(series::TruncatedBSeries)

The maximal `order` of a rooted tree with non-vanishing coefficient in the
truncated B-series `series`.

See also [`order_of_accuracy`](@ref).
"""
RootedTrees.order(series::TruncatedBSeries) = order(series.coef.keys[end])

# vector space interface
for op in (:+, :-)
    @eval function ($op)(series1::TruncatedBSeries{T},
                         series2::TruncatedBSeries{T}) where {T}
        # truncate the longer B-series
        if order(series1) > order(series2)
            series_keys = keys(series1)
        else
            series_keys = keys(series2)
        end

        V = promote_type(valtype(series1), valtype(series2))
        series_result = TruncatedBSeries{T, V}()
        for (key, val1, val2) in zip(series_keys, values(series1), values(series2))
            series_result[key] = ($op)(val1, val2)
        end

        return series_result
    end

    @eval function ($op)(series::TruncatedBSeries)
        series_keys = keys(series)

        T = keytype(series)
        V = valtype(series)
        series_result = TruncatedBSeries{T, V}()
        for (key, val) in zip(series_keys, values(series))
            series_result[key] = ($op)(val)
        end

        return series_result
    end
end

for op in (:*, :/)
    @eval function ($op)(series::TruncatedBSeries, scalar::Number)
        series_keys = keys(series)

        T = keytype(series)
        V = promote_type(valtype(series), typeof(scalar))
        series_result = TruncatedBSeries{T, V}()
        for (key, val) in zip(series_keys, values(series))
            series_result[key] = ($op)(val, scalar)
        end

        return series_result
    end
end

for op in (:*, :\)
    @eval function ($op)(scalar::Number, series::TruncatedBSeries)
        series_keys = keys(series)

        T = keytype(series)
        V = promote_type(valtype(series), typeof(scalar))
        series_result = TruncatedBSeries{T, V}()
        for (key, val) in zip(series_keys, values(series))
            series_result[key] = ($op)(scalar, val)
        end

        return series_result
    end
end

# Return a function returning an iterator over all rooted trees used as
# keys for the B-series when given an order of the trees.
function _iterator_type(::TruncatedBSeries{<:RootedTree})
    return RootedTreeIterator
end

function _iterator_type(::TruncatedBSeries{<:BicoloredRootedTree})
    return BicoloredRootedTreeIterator
end

"""
    ExactSolution{V}()

Lazy representation of the B-series of the exact solution of an ordinary
differential equation using coefficients of type at least as representative as
`V`.
"""
struct ExactSolution{V} end

function Base.getindex(::ExactSolution{V}, t::AbstractRootedTree) where {V}
    convert(V, 1 // 1) / γ(t)
end

# general interface methods of iterators for `ExactSolution`
Base.IteratorSize(::Type{<:ExactSolution}) = Base.SizeUnknown()
Base.eltype(::Type{ExactSolution{V}}) where {V} = V
Base.valtype(exact::ExactSolution) = valtype(typeof(exact))
Base.valtype(::Type{ExactSolution{V}}) where {V} = V

function Base.iterate(exact::ExactSolution)
    iterator = RootedTreeIterator(0)
    t, state = iterate(iterator)
    (exact[t], (state, iterator))
end

function Base.iterate(exact::ExactSolution, state_iterator)
    state, iterator = state_iterator
    t_state = iterate(iterator, state)
    if t_state === nothing
        iterator = RootedTreeIterator(iterator.order + 1)
        t_state = iterate(iterator)
    end
    t, state = t_state
    (exact[t], (state, iterator))
end

# internal interface of B-series
# @inline evaluation_type(::ExactSolution) = LazyEvaluation() # this is the default assumption

"""
    ExactSolution(series_integrator)

A representation of the B-series of the exact solution of an ODE using the same
type of coefficients as the B-series `series_integrator`.
"""
function ExactSolution(series_integrator)
    ExactSolution(series_integrator, evaluation_type(series_integrator))
end

function ExactSolution(series_integrator, ::LazyEvaluation)
    ExactSolution{valtype(series_integrator)}()
end

function ExactSolution(series_integrator, ::EagerEvaluation)
    exact = ExactSolution{valtype(series_integrator)}()
    T = keytype(series_integrator)
    V = valtype(series_integrator)
    series = TruncatedBSeries{T, V}()
    series_keys = keys(series_integrator)

    iter = iterate(series_keys)
    if iter !== nothing
        t, t_state = iter
        if isempty(t)
            series[t] = one(V)
            iter = iterate(series_keys, t_state)
        end
    end

    while iter !== nothing
        t, t_state = iter
        series[t] = exact[t]
        iter = iterate(series_keys, t_state)
    end
    series
end

# vector space interface
for op in (:+, :-)
    @eval function ($op)(series1::TruncatedBSeries, series2::ExactSolution)
        # truncate the B-series of the exact solution
        series_keys = keys(series1)

        T = keytype(series1)
        V = promote_type(valtype(series1), valtype(series2))
        series_result = TruncatedBSeries{T, V}()
        for (key, val1, val2) in zip(series_keys, values(series1), values(series2))
            series_result[key] = ($op)(val1, val2)
        end

        return series_result
    end

    @eval function ($op)(series1::ExactSolution, series2::TruncatedBSeries)
        # truncate the B-series of the exact solution
        series_keys = keys(series2)

        T = keytype(series2)
        V = promote_type(valtype(series1), valtype(series2))
        series_result = TruncatedBSeries{T, V}()
        for (key, val1, val2) in zip(series_keys, values(series1), values(series2))
            series_result[key] = ($op)(val1, val2)
        end

        return series_result
    end
end

# investigate properties of B-series
"""
    order_of_accuracy(series; kwargs...)

Determine the order of accuracy of the B-series `series`. By default, the
comparison with the coefficients of the exact solution is performed using
`isequal`. If keyword arguments such as absolute/relative tolerances `atol`/`rtol`
are given or floating point numbers are used, the comparison is performed using
`isapprox` and the keyword arguments `kwargs...` are forwarded.

See also [`order`](@ref), [`ExactSolution`](@ref).
"""
function order_of_accuracy(series::TruncatedBSeries; kwargs...)
    if isempty(kwargs) && !(valtype(series) <: AbstractFloat)
        compare = isequal
    else
        compare = (a, b) -> isapprox(a, b; kwargs...)
    end

    exact = ExactSolution{valtype(series)}()
    iterator = _iterator_type(series)

    for o in 0:order(series)
        # Iterate over all rooted trees used as keys in `series`
        # of a given order `o`.
        for t in iterator(o)
            if compare(series[t], exact[t]) == false
                return order(t) - 1
            end
        end
    end

    return order(series)
end

# construct B-series
"""
    bseries(f::Function, order, iterator_type=RootedTreeIterator)

Return a truncated B-series up to the specified `order` with coefficients
determined by `f`. The type of rooted trees is determined by the `iterator_type`,
which can be `RootedTreeIterator` or `BicoloredRootedTreeIterator`. Calling
`f(t, series)` needs to return the coefficient of the rooted tree `t` of the
desired series in a type-stable manner. For the empty tree, `f` is called as
`f(t, nothing)`. Otherwise, the series constructed so far is passed as second
argument, allowing one to access values of lower-order trees.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be
    multiplied by a power of the time step divided by the `symmetry` of the
    rooted tree and multiplied by the corresponding elementary differential
    of the input vector field ``f``.
    See also [`evaluate`](@ref).

# Examples

The B-series of the average vector field (AVF) method is given by
``b(.) = 1`` and ``b([t_1, ..., t_n]) = b(t_1)...b(t_n) / (n + 1)``, see
- Elena Celledoni, Robert I. McLachlan, David I. McLaren, Brynjulf Owren,
  G. Reinout W. Quispel, and William M. Wright.
  "Energy-preserving Runge-Kutta methods."
  ESAIM: Mathematical Modelling and Numerical Analysis 43, no. 4 (2009): 645-649.
  [DOI: 10.1051/m2an/2009020](https://doi.org/10.1051/m2an/2009020)

We can generate this as follows.
```jldoctest
julia> series = bseries(3) do t, series
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
TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Rational{Int64}} with 5 entries:
  RootedTree{Int64}: Int64[]   => 1//1
  RootedTree{Int64}: [1]       => 1//1
  RootedTree{Int64}: [1, 2]    => 1//2
  RootedTree{Int64}: [1, 2, 3] => 1//4
  RootedTree{Int64}: [1, 2, 2] => 1//3
```
"""
function bseries(f::Function, order, iterator_type = RootedTreeIterator)
    # Get the coefficient of the empty tree
    t = first(iterator_type(0))
    v = f(t, nothing)

    V_tmp = typeof(v)
    if V_tmp <: Integer
        # If people use integer coefficients, they will likely want to have results
        # as exact as possible. However, general terms are not integers. Thus, we
        # use rationals instead.
        V = Rational{V_tmp}
    else
        V = V_tmp
    end

    # Setup the series
    series = TruncatedBSeries{typeof(t), V}()
    series[copy(t)] = v

    for o in 1:order
        for t in iterator_type(o)
            series[copy(t)] = f(t, series)
        end
    end

    return series
end

"""
    bseries(rk::RungeKuttaMethod, order)
    bseries(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)

Compute the B-series of the Runge-Kutta method `rk` with Butcher coefficients
`A, b, c` up to a prescribed integer `order`.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be
    multiplied by a power of the time step divided by the `symmetry` of the
    rooted tree and multiplied by the corresponding elementary differential
    of the input vector field ``f``.
    See also [`evaluate`](@ref).
"""
function bseries(rk::RungeKuttaMethod, order)
    V_tmp = eltype(rk)
    if V_tmp <: Integer
        # If people use integer coefficients, they will likely want to have results
        # as exact as possible. However, general terms are not integers. Thus, we
        # use rationals instead.
        V = Rational{V_tmp}
    else
        V = V_tmp
    end
    series = TruncatedBSeries{RootedTree{Int, Vector{Int}}, V}()

    series[rootedtree(Int[])] = one(V)
    for o in 1:order
        for t in RootedTreeIterator(o)
            series[copy(t)] = elementary_weight(t, rk)
        end
    end

    return series
end

function bseries(A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                 order)
    rk = RungeKuttaMethod(A, b, c)
    bseries(rk, order)
end

# TODO: bseries(rk::RungeKuttaMethod)
# should create a lazy version, optionally a memoized one

"""
    AverageVectorFieldMethod([T=Rational{Int}])

Construct a representation of the average vector field (AVF) method using
coefficients of type `T`. You can pass it as argument to [`bseries`](@ref)
to construct the corresponding B-series.

# Examples

We can generate this as follows.
```jldoctest
julia> series = bseries(AverageVectorFieldMethod(), 3)
TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Rational{Int64}} with 5 entries:
  RootedTree{Int64}: Int64[]   => 1//1
  RootedTree{Int64}: [1]       => 1//1
  RootedTree{Int64}: [1, 2]    => 1//2
  RootedTree{Int64}: [1, 2, 3] => 1//4
  RootedTree{Int64}: [1, 2, 2] => 1//3
```

# References

The B-series of the average vector field (AVF) method is given by
``b(.) = 1`` and ``b([t_1, ..., t_n]) = b(t_1)...b(t_n) / (n + 1)``, see

- Elena Celledoni, Robert I. McLachlan, David I. McLaren, Brynjulf Owren,
  G. Reinout W. Quispel, and William M. Wright.
  "Energy-preserving Runge-Kutta methods."
  ESAIM: Mathematical Modelling and Numerical Analysis 43, no. 4 (2009): 645-649.
  [DOI: 10.1051/m2an/2009020](https://doi.org/10.1051/m2an/2009020)
"""
struct AverageVectorFieldMethod{T} end

AverageVectorFieldMethod() = AverageVectorFieldMethod(Rational{Int})
AverageVectorFieldMethod(::Type{T}) where {T} = AverageVectorFieldMethod{T}()

"""
    bseries(avf::AverageVectorFieldMethod, order)

Compute the B-series of the [`AverageVectorFieldMethod`](@ref) up to a
prescribed integer `order`.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by `bseries` need to be
    multiplied by a power of the time step divided by the `symmetry` of the
    rooted tree and multiplied by the corresponding elementary differential
    of the input vector field ``f``.
    See also [`evaluate`](@ref).
"""
function bseries(::AverageVectorFieldMethod{T}, o) where {T}
    bseries(o) do t, series
        if order(t) in (0, 1)
            return one(T)
        else
            v = one(T)
            n = 0
            for subtree in SubtreeIterator(t)
                v *= series[subtree]
                n += 1
            end
            return v / (n + 1)
        end
    end
end

"""
    bseries(ark::AdditiveRungeKuttaMethod, order)

Compute the B-series of the additive Runge-Kutta method `ark` up to a prescribed
integer `order`.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be
    multiplied by a power of the time step divided by the `symmetry` of the
    colored rooted tree and multiplied by the corresponding elementary
    differential of the input vector fields ``f^\\nu``.
    See also [`evaluate`](@ref).
"""
function bseries(ark::AdditiveRungeKuttaMethod, order)
    if length(ark.rks) != 2
        throw(ArgumentError("Only AdditiveRungeKuttaMethod with a dual splitting are supported. Got an ARK with $(length(ark.rks)) methods."))
    end

    V_tmp = eltype(ark)
    if V_tmp <: Integer
        # If people use integer coefficients, they will likely want to have results
        # as exact as possible. However, general terms are not integers. Thus, we
        # use rationals instead.
        V = Rational{V_tmp}
    else
        V = V_tmp
    end
    series = TruncatedBSeries{BicoloredRootedTree{Int, Vector{Int}, Vector{Bool}}, V}()

    series[rootedtree(Int[], Bool[])] = one(V)
    for o in 1:order
        for t in BicoloredRootedTreeIterator(o)
            series[copy(t)] = elementary_weight(t, ark)
        end
    end

    return series
end

# TODO: bseries(ark::AdditiveRungeKuttaMethod)
# should create a lazy version, optionally a memoized one

"""
    bseries(ros::RosenbrockMethod, order)

Compute the B-series of the Rosenbrock method `ros` up to a prescribed
integer `order`.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be
    multiplied by a power of the time step divided by the `symmetry` of the
    rooted tree and multiplied by the corresponding elementary differential
    of the input vector field ``f``.
    See also [`evaluate`](@ref).
"""
function bseries(ros::RosenbrockMethod, order)
    V_tmp = eltype(ros)
    if V_tmp <: Integer
        # If people use integer coefficients, they will likely want to have results
        # as exact as possible. However, general terms are not integers. Thus, we
        # use rationals instead.
        V = Rational{V_tmp}
    else
        V = V_tmp
    end
    series = TruncatedBSeries{RootedTree{Int, Vector{Int}}, V}()

    series[rootedtree(Int[])] = one(V)
    for o in 1:order
        for t in RootedTreeIterator(o)
            series[copy(t)] = elementary_weight(t, ros)
        end
    end

    return series
end

# TODO: bseries(ros::RosenbrockMethod)
# should create a lazy version, optionally a memoized one

"""
    ContinuousStageRungeKuttaMethod(M)

A struct that describes a continuous stage Runge-Kutta (CSRK) method. It can
be constructed by passing the parameter matrix `M` as described by Miyatake and
Butcher (2016). You can compute the B-series of the method by [`bseries`](@ref).

# References

- Yuto Miyatake and John C. Butcher.
  "A characterization of energy-preserving methods and the construction of
  parallel integrators for Hamiltonian systems."
  SIAM Journal on Numerical Analysis 54, no. 3 (2016):
  [DOI: 10.1137/15M1020861](https://doi.org/10.1137/15M1020861)
"""
struct ContinuousStageRungeKuttaMethod{MatT <: AbstractMatrix}
    matrix::MatT
end

# TODO: Documentation (tutorial), Base.show, etc.

"""
    bseries(csrk::ContinuousStageRungeKuttaMethod, order)

Compute the B-series of the [`ContinuousStageRungeKuttaMethod`](@ref) `csrk`
up to the prescribed integer `order` as described by Miyatake & Butcher (2016).

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be
    multiplied by a power of the time step divided by the `symmetry` of the
    rooted tree and multiplied by the corresponding elementary differential
    of the input vector field ``f``.
    See also [`evaluate`](@ref).

# Examples

The [`AverageVectorFieldMethod`](@ref) is given by the parameter matrix with
single entry one.

```jldoctest
julia> M = fill(1//1, 1, 1)
1×1 Matrix{Rational{Int64}}:
 1//1

julia> series = bseries(ContinuousStageRungeKuttaMethod(M), 4)
TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Rational{Int64}} with 9 entries:
  RootedTree{Int64}: Int64[]      => 1//1
  RootedTree{Int64}: [1]          => 1//1
  RootedTree{Int64}: [1, 2]       => 1//2
  RootedTree{Int64}: [1, 2, 3]    => 1//4
  RootedTree{Int64}: [1, 2, 2]    => 1//3
  RootedTree{Int64}: [1, 2, 3, 4] => 1//8
  RootedTree{Int64}: [1, 2, 3, 3] => 1//6
  RootedTree{Int64}: [1, 2, 3, 2] => 1//6
  RootedTree{Int64}: [1, 2, 2, 2] => 1//4

julia> series - bseries(AverageVectorFieldMethod(), order(series))
TruncatedBSeries{RootedTree{Int64, Vector{Int64}}, Rational{Int64}} with 9 entries:
  RootedTree{Int64}: Int64[]      => 0//1
  RootedTree{Int64}: [1]          => 0//1
  RootedTree{Int64}: [1, 2]       => 0//1
  RootedTree{Int64}: [1, 2, 3]    => 0//1
  RootedTree{Int64}: [1, 2, 2]    => 0//1
  RootedTree{Int64}: [1, 2, 3, 4] => 0//1
  RootedTree{Int64}: [1, 2, 3, 3] => 0//1
  RootedTree{Int64}: [1, 2, 3, 2] => 0//1
  RootedTree{Int64}: [1, 2, 2, 2] => 0//1
```

# References

- Yuto Miyatake and John C. Butcher.
  "A characterization of energy-preserving methods and the construction of
  parallel integrators for Hamiltonian systems."
  SIAM Journal on Numerical Analysis 54, no. 3 (2016):
  [DOI: 10.1137/15M1020861](https://doi.org/10.1137/15M1020861)
"""
function bseries(csrk::ContinuousStageRungeKuttaMethod, order)
    bseries(order) do t, series
        elementary_weight(t, csrk)
    end
end

# TODO: bseries(csrk::ContinuousStageRungeKuttaMethod)
# should create a lazy version, optionally a memoized one

"""
    elementary_weight(t::RootedTree, csrk::ContinuousStageRungeKuttaMethod)

Compute the elementary weight Φ(`t`) of the
[`ContinuousStageRungeKuttaMethod`](@ref) `csrk`
for a rooted tree `t``.

# References

- Yuto Miyatake and John C. Butcher.
  "A characterization of energy-preserving methods and the construction of
  parallel integrators for Hamiltonian systems."
  SIAM Journal on Numerical Analysis 54, no. 3 (2016):
  [DOI: 10.1137/15M1020861](https://doi.org/10.1137/15M1020861)
"""
function RootedTrees.elementary_weight(t::RootedTree,
                                       csrk::ContinuousStageRungeKuttaMethod)
    # See Miyatake & Butcher (2016), p. 1998, right before eq. (2.8)
    if order(t) <= 1
        return one(eltype(csrk.matrix))
    end

    # First, we compute the coefficient matrix `A_τζ` of the method from
    # the matrix `M = csrk.matrix`. We compute it only once and pass it to
    # `derivative_weight` below to save additional computations.
    A_τζ = _matrix_a(csrk)

    # The elementary weight Φ is given as
    #   Φ(t) = ∫₀¹ B_τ ϕ_τ(t₁) ... ϕ_τ(tₘ) dτ
    # where
    #   B_ζ = A_1ζ
    # Since the polynomial matrix `A_τζ` is stored as a polyomial in ζ
    # with coefficients as polynomials in τ, setting `τ = 1` means that
    # we need to evaluate all coefficients at 1 and interpret the result
    # as a polynomial in τ.
    integrand = let
        v = map(p -> p(1), Polynomials.coeffs(A_τζ))
        Polynomial{eltype(v), :τ}(v)
    end
    # Now, we can loop over all subtrees and simply update the integrand.
    for subtree in SubtreeIterator(t)
        ϕ = derivative_weight(subtree, A_τζ, csrk)
        integrand = integrand * ϕ
    end

    # The antiderivative comes with a zero constant term. Thus, we can just
    # evaluate it at 1 to get the value of the integral from 0 to 1.
    antiderivative = Polynomials.integrate(integrand)
    result = antiderivative(1)

    return result
end

# See Miyatake & Butcher (2016), p. 1998, right before eq. (2.8)
function derivative_weight(t::RootedTree, A_τζ, csrk::ContinuousStageRungeKuttaMethod)

    # The derivative weight ϕ_τ is given as
    #   ϕ_τ(t) = ∫₀¹ A_τζ ϕ_ζ(t₁) ... ϕ_ζ(tₘ) dζ
    # Since the polynomial matrix `A_τζ` is stored as a polyomial in ζ
    # with coefficients as polynomials in τ, we need to interpret the
    # return value of the inner `derivative_weight` as the constant
    # polynomial in ζ with coefficients given as polynomials in τ.
    integrand = A_τζ
    for subtree in SubtreeIterator(t)
        ϕ = derivative_weight(subtree, A_τζ, csrk)
        integrand = integrand * Polynomial{typeof(ϕ), :ζ}(Polynomials.coeffs(ϕ))
    end

    # The antiderivative comes with a zero constant term. Thus, we can just
    # evaluate it at 1 to get the value of the integral from 0 to 1.
    antiderivative = Polynomials.integrate(integrand)
    result = antiderivative(1)

    return result
end

# Equation (2.6) of Miyatake & Butcher (2016)
function _matrix_a(csrk::ContinuousStageRungeKuttaMethod)
    M = csrk.matrix
    s = size(M, 1)
    T_tmp = eltype(M)
    if T_tmp <: Integer
        # If people use integer coefficients, they will likely want to have results
        # as exact as possible. However, general terms are not integers. Thus, we
        # use rationals instead.
        T = Rational{T_tmp}
    else
        T = T_tmp
    end

    τ = Vector{Polynomial{T, :τ}}(undef, s)
    for i in 1:s
        v = zeros(T, i + 1)
        v[end] = 1 // i
        τ[i] = Polynomial{T, :τ}(v)
    end

    Mτ = M' * τ
    A_τζ = Polynomial{eltype(Mτ), :ζ}(Mτ)

    return A_τζ
end

# TODO: Documentation (tutorial), Base.show, etc.
"""
    MultirateInfinitesimalSplitMethod(A, D, G, c)

# References

- Knoth, Oswald, and Joerg Wensch.
  "Generalized split-explicit Runge-Kutta methods for the
  compressible Euler equations".
  Monthly Weather Review 142, no. 5 (2014): 2067-2081.
  [DOI: 10.1175/MWR-D-13-00068.1](https://doi.org/10.1175/MWR-D-13-00068.1)

!!! warning "Experimental code"
    This code is considered to be experimental at the moment
    and can change any time.
"""
struct MultirateInfinitesimalSplitMethod{T, PolyMatT <: AbstractMatrix{<:Polynomial{T}},
                                         MatT <: AbstractMatrix{T},
                                         VecT <: AbstractVector{T}} <:
       RootedTrees.AbstractTimeIntegrationMethod
    A::PolyMatT
    D::MatT
    G::MatT
    c::VecT
end

# TODO: Deduce `c` from other parameters?
function MultirateInfinitesimalSplitMethod(A::AbstractMatrix{<:Polynomial},
                                           D::AbstractMatrix,
                                           G::AbstractMatrix,
                                           c::AbstractVector)
    T = promote_type(eltype(eltype(A)), eltype(D), eltype(G), eltype(c))
    PolyT = typeof(zero(first(A)) + zero(T))
    _A = PolyT.(A)
    _D = T.(D)
    _G = T.(G)
    _c = T.(c)
    return MultirateInfinitesimalSplitMethod(_A, _D, _G, _c)
end

Base.eltype(mis::MultirateInfinitesimalSplitMethod{T}) where {T} = T

"""
    bseries(mis::MultirateInfinitesimalSplitMethod, order)

Compute the B-series of the multirate infinitesimal split method `mis`
up to a prescribed integer `order`.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be
    multiplied by a power of the time step divided by the `symmetry` of the
    colored rooted tree and multiplied by the corresponding elementary
    differential of the input vector fields ``f^\\nu``.
    See also [`evaluate`](@ref).
"""
function bseries(mis::MultirateInfinitesimalSplitMethod, order)
    V_tmp = eltype(mis)
    if V_tmp <: Integer
        # If people use integer coefficients, they will likely want to have results
        # as exact as possible. However, general terms are not integers. Thus, we
        # use rationals instead.
        V = Rational{V_tmp}
    else
        V = V_tmp
    end

    prototype_scalar = TruncatedBSeries{BicoloredRootedTree{Int, Vector{Int}, Vector{Bool}},
                                        V}()
    prototype_scalar[rootedtree(Int[], Bool[])] = one(V)

    poly_one = Polynomial{V, :x}([one(V)])
    poly_zero = zero(poly_one)
    prototype_polynomial = TruncatedBSeries{
                                            BicoloredRootedTree{Int, Vector{Int},
                                                                Vector{Bool}},
                                            typeof(poly_one)}()
    prototype_polynomial[rootedtree(Int[], Bool[])] = poly_one

    A = mis.A
    D = mis.D
    G = mis.G
    ns = size(A, 2)
    Z = Vector{typeof(prototype_polynomial)}(undef, ns + 1)
    η = Vector{typeof(prototype_scalar)}(undef, ns + 1)
    for i in 1:(ns + 1)
        Z[i] = copy(prototype_polynomial)
        η[i] = copy(prototype_scalar)
    end

    for o in 1:order
        for t_ in BicoloredRootedTreeIterator(o)
            t = copy(t_)
            for i in 1:(ns + 1)
                phi = poly_zero
                for j in 1:(i - 1)
                    r = poly_one
                    for subtree in RootedTrees.SubtreeIterator(t)
                        if subtree.color_sequence[1]
                            r = r * Z[i][subtree]
                        else
                            r = r * η[j][subtree]
                        end
                    end

                    v = [one(V)]
                    for k in 0:(length(A[i, j]) - 1)
                        phi = phi + Polynomial{V, :x}(v) * A[i, j][k] * r
                        v = [0; v]
                    end
                end

                for j in 1:(i - 1)
                    phi = phi + G[i, j] * (η[j][t] - η[1][t])
                end

                phi = Polynomials.integrate(phi)
                phi[0] = 0
                for j in 1:(i - 1)
                    phi[0] = phi[0] + D[i, j] * (η[j][t] - η[1][t])
                end

                Z[i][t] = phi
                η[i][t] = phi(1)
            end
        end
    end

    series = η[end]
    return series
end

# TODO: bseries(mis::MultirateInfinitesimalSplitMethod)
# should create a lazy version, optionally a memoized one

"""
    substitute(b, a, t::AbstractRootedTree)

Compute the coefficient corresponding to the tree `t` of the B-series that is
formed by substituting the B-series `b` into the B-series `a`. It is assumed
that the B-series `b` has the coefficient zero of the empty tree.

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function substitute(b, a, t::AbstractRootedTree)
    result = zero(first(values(a)) * first(values(b)))

    for (forest, skeleton) in PartitionIterator(t)
        update = a[skeleton]
        update isa Rational && iszero(update) && continue
        for tree in forest
            update *= b[tree]
        end
        result += update
    end

    return result
end

"""
    substitute(b, a)

Substitute the B-series `b` into the B-series `a`. It is assumed that the
B-series `b` has the coefficient zero of the empty tree.

In the notation of Chartier, Hairer and Vilmart (2010), we have
`substitute(b, a) = b ★ a`.

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function substitute(b, a)
    series_keys = keys(b)
    series = empty(b)

    t = first(series_keys)
    @assert isempty(t)
    series[t] = a[t]

    for t in Iterators.drop(series_keys, 1)
        coefficient = substitute(b, a, t)
        series[t] = coefficient
    end

    return series
end

"""
    compose(b, a, t::RootedTree)

Compute the coefficient corresponding to the tree `t` of the B-series that is
formed by composing the B-series `a` with the B-series `b`. It is assumed that
the B-series `b` has the coefficient unity of the empty tree.

# References

Section 3.1 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function compose(b, a, t::RootedTree)
    result = zero(first(values(a)) * first(values(b)))

    for (forest, subtree) in SplittingIterator(t)
        update = a[subtree]
        update isa Rational && iszero(update) && continue
        for tree in forest
            update *= b[tree]
        end
        result += update
    end

    return result
end

"""
    compose(b, a; normalize_stepsize=false)

Compose the B-series `a` with the B-series `b`. It is assumed that the B-series
`b` has the coefficient unity of the empty tree.

In the notation of Chartier, Hairer and Vilmart (2010), we have
`compose(b, a) = b ⋅ a`. Note that this means that method `b` is applied first,
followed by method `a`.

If `normalize_stepsize = true`, the coefficients of the returned B-series
are divided by `2^order(t)` for each rooted tree `t`. This normalizes the step
size so that the resulting numerical integrator B-series uses the same step size
as the input series (instead of a doubled step size).

# References

Section 3.1 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function compose(b, a; normalize_stepsize = false)
    series_keys = keys(b)
    series = empty(b)

    for t in series_keys
        coefficient = compose(b, a, t)
        if normalize_stepsize
            coefficient /= 2^order(t)
        end
        series[t] = coefficient
    end

    return series
end

"""
    compose(b1, b2, bs...; normalize_stepsize=false)

Compose the B-series `b1`, `b2`, `bs...`. It is assumed that all B-series
have the coefficient unity of the empty tree.

In the notation of Chartier, Hairer and Vilmart (2010), we have
`compose(b1, b2, b3) = b1 ⋅ b2 ⋅ b3`. Note that this product is associative and
has to be read from left to right, i.e., method `b1` is applied first, followed
by `b2`, `bs...`.

If `normalize_stepsize = true`, the coefficients of the returned B-series
are divided by `n^order(t)` for each rooted tree `t`, where `n` is the total
number of composed B-series. This normalizes the step size so that the resulting
numerical integrator B-series uses the same step size as the input series
(instead of an `n`-fold step size).

# References

Section 3.1 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function compose(b1, b2, bs::Vararg{Any, N}; normalize_stepsize = false) where {N}
    series = compose(b1, b2; normalize_stepsize = false)

    for b in bs
        series = compose(series, b; normalize_stepsize = false)
    end

    if normalize_stepsize
        for (t, v) in series
            series[t] = v / (N + 2)^order(t)
        end
    end

    return series
end

"""
    evaluate(f, u, dt, series, reduce_order_by=0)

Evaluate the B-series `series` specialized to the ordinary differential equation
``u'(t) = f(u(t))`` with vector field `f` and dependent variables `u` for a
time step size `dt`.

Here, `u` is assumed to be a vector of symbolic variables and `f` is assumed
to be a vector of expressions in these variables for plain B-series. For
B-series with colored trees, `f` must be a tuple of vectors of expressions in
the variables `u`. Currently, symbolic variables from

- [SymEngine.jl](https://github.com/symengine/SymEngine.jl),
- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), and
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

are supported.

The powers of `dt` can be controlled by `reduce_order_by` to make them different
from the usual `order(t)` for a rooted tree `t`. This can be useful in the
context of [`modified_equation`](@ref)s or [`modifying_integrator`](@ref)s,
where the B-series coefficients are those of ``h fₕ``, i.e., they contain an
additional power of `dt`. In this case, the B-series of the vector field can
be obtained using `reduce_order_by = 1`.

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function evaluate(f, u, dt, series, reduce_order_by = 0)
    _evaluate(f, u, dt, series, evaluation_type(series), reduce_order_by)
end

function _evaluate(f, u, dt, series, ::EagerEvaluation, reduce_order_by)
    differentials = elementary_differentials(f, u, order(series))

    # An additive decomposition is indicated by a tuple of vectors. A single
    # vector field is assumed to be a vector, as everywhere else.
    if f isa NTuple{N, AbstractVector} where {N}
        result = zero(first(f))
    else
        result = zero(f)
    end

    for t in keys(series)
        # Otherwise, SymPy.jl might result in
        #   DomainError with -1:
        #   Cannot raise an integer x to a negative power -1.
        #   Convert input to float.
        power = order(t) - reduce_order_by
        if power > 0
            dt_term = dt^power
        elseif power == 0
            dt_term = one(dt)
        else
            dt_term = inv(dt^(-power))
        end

        result += dt_term / symmetry(t) * series[t] * differentials[t]
    end
    result
end

"""
    modified_equation(series_integrator)

Compute the B-series of the modified equation of the time integration method
with B-series `series_integrator`.

Given an ordinary differential equation (ODE) ``u'(t) = f(u(t))`` and a
Runge-Kutta method, the idea is to interpret the numerical solution with
given time step size as exact solution of a modified ODE ``u'(t) = fₕ(u(t))``.
This method returns the B-series of ``h fₕ``.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be
    multiplied by a power of the time step divided by the `symmetry` of the
    rooted tree and multiplied by the corresponding elementary differential
    of the input vector field ``f``.
    See also [`evaluate`](@ref).

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function modified_equation(series_integrator)
    _modified_equation(series_integrator, evaluation_type(series_integrator))
end

function _modified_equation(series_integrator, ::EagerEvaluation)
    V = valtype(series_integrator)

    # B-series of the exact solution
    # We could just use the lazy version
    #   series_ex = ExactSolution{V}()
    # However, we need to access elements of `series_ex` more than once in the
    # substitution below. Thus, it's cheaper to compute every entry only once and
    # re-use it later.
    series_ex = ExactSolution(series_integrator)

    # Prepare B-series of the modified equation
    series_keys = keys(series_integrator)
    series = empty(series_integrator)
    for t in series_keys
        series[t] = zero(V)
    end

    iter = iterate(series_keys)
    if iter !== nothing
        t, t_state = iter
        if isempty(t)
            iter = iterate(series_keys, t_state)
            if iter !== nothing
                t, t_state = iter
            end
        end

        series[t] = series_integrator[t]
        iter = iterate(series_keys, t_state)
    end

    # Recursively solve
    #   substitute(series, series_ex, t) == series_integrator[t]
    # This works because
    #   substitute(series, series_ex, t) = series[t] + lower order terms
    # Since the `keys` are ordered, we don't need to use nested loops of the form
    #   for o in 2:order
    #     for _t in RootedTreeIterator(o)
    #       t = copy(_t)
    # which are slightly less efficient due to additional computations and
    # allocations.
    while iter !== nothing
        t, t_state = iter
        series[t] += series_integrator[t] - substitute(series, series_ex, t)
        iter = iterate(series_keys, t_state)
    end

    return series
end

"""
    modified_equation(rk::RungeKuttaMethod, order)
    modified_equation(A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                      order)

Compute the B-series of the [`modified_equation`](@ref) of the Runge-Kutta
method `rk` with Butcher coefficients `A, b, c` up to the prescribed `order`.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be
    multiplied by a power of the time step divided by the `symmetry` of the
    rooted tree and multiplied by the corresponding elementary differential
    of the input vector field ``f``.
    See also [`evaluate`](@ref).
"""
function modified_equation(rk::RungeKuttaMethod, order)
    # B-series of the Runge-Kutta method
    series = bseries(rk, order)
    modified_equation(series)
end

function modified_equation(A::AbstractMatrix, b::AbstractVector,
                           c::AbstractVector, order)
    rk = RungeKuttaMethod(A, b, c)
    modified_equation(rk, order)
end

"""
    modified_equation(f, u, dt, series_integrator)

Compute the B-series of the [`modified_equation`](@ref) of the time integration
method with B-series `series_integrator` with respect to the ordinary
differential equation ``u'(t) = f(u(t))`` with vector field `f` and dependent
variables `u` for a time step size `dt`.

Here, `u` is assumed to be a vector of symbolic variables and `f` is assumed
to be a vector of expressions in these variables for plain B-series. For
B-series with colored trees, `f` must be a tuple of vectors of expressions in
the variables `u`. Currently, symbolic variables from

- [SymEngine.jl](https://github.com/symengine/SymEngine.jl),
- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), and
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

are supported.
"""
function modified_equation(f, u, dt, series_integrator)
    series = modified_equation(series_integrator)
    reduce_order_by = 1 # reduce the powers of `dt` by one since dt*f is given
    # by recursively solving `substitute`
    evaluate(f, u, dt, series, reduce_order_by)
end

"""
    modified_equation(f, u, dt, rk::RungeKuttaMethod, order)
    modified_equation(f, u, dt,
                      A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                      order)

Compute the B-series of the [`modified_equation`](@ref) of the Runge-Kutta
method `rk` with Butcher coefficients `A, b, c` up to the prescribed `order` with
respect to the ordinary differential equation ``u'(t) = f(u(t))`` with vector
field `f` and dependent variables `u` for a time step size `dt`.

Here, `u` is assumed to be a vector of symbolic variables and `f` is assumed
to be a vector of expressions in these variables for plain B-series. For
B-series with colored trees, `f` must be a tuple of vectors of expressions in
the variables `u`. Currently, symbolic variables from

- [SymEngine.jl](https://github.com/symengine/SymEngine.jl),
- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), and
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

are supported.
"""
function modified_equation(f, u, dt, rk::RungeKuttaMethod, order)
    series_integrator = bseries(rk, order)
    modified_equation(f, u, dt, series_integrator)
end

function modified_equation(f, u, dt,
                           A::AbstractMatrix, b::AbstractVector,
                           c::AbstractVector, order)
    rk = RungeKuttaMethod(A, b, c)
    modified_equation(f, u, dt, rk, order)
end

"""
    modifying_integrator(series_integrator)

Compute the B-series of a "modifying integrator" equation of the time
integration method with B-series `series_integrator`.

Given an ordinary differential equation (ODE) ``u'(t) = f(u(t))`` and a
Runge-Kutta method, the idea is to find a modified ODE ``u'(t) = fₕ(u(t))``
such that the numerical solution with given time step size is the exact solution
of the original ODE. This method returns the B-series of ``h fₕ``.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be
    multiplied by a power of the time step divided by the `symmetry` of the
    rooted tree and multiplied by the corresponding elementary differential
    of the input vector field ``f``.
    See also [`evaluate`](@ref).

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function modifying_integrator(series_integrator)
    _modifying_integrator(series_integrator, evaluation_type(series_integrator))
end

function _modifying_integrator(series_integrator, ::EagerEvaluation)
    V = valtype(series_integrator)

    # B-series of the exact solution
    # Since we access each coefficient of this B-series only once below, it's
    # okay to use the lazy version of the exact solution.
    series_ex = ExactSolution{V}()

    # Prepare B-series of the modifying integrator equation
    series_keys = keys(series_integrator)
    series = empty(series_integrator)
    for t in series_keys
        series[t] = zero(V)
    end

    iter = iterate(series_keys)
    if iter !== nothing
        t, t_state = iter
        if isempty(t)
            iter = iterate(series_keys, t_state)
            if iter !== nothing
                t, t_state = iter
            end
        end

        series[t] = series_integrator[t]
        iter = iterate(series_keys, t_state)
    end

    # Recursively solve
    #   substitute(series, series_integrator, t) == series_ex[t]
    # This works because
    #   substitute(series, series_integrator, t) = series[t] + lower order terms
    # Since the `keys` are ordered, we don't need to use nested loops of the form
    #   for o in 2:order
    #     for _t in RootedTreeIterator(o)
    #       t = copy(_t)
    # which are slightly less efficient due to additional computations and
    # allocations.
    while iter !== nothing
        t, t_state = iter
        series[t] += series_ex[t] - substitute(series, series_integrator, t)
        iter = iterate(series_keys, t_state)
    end

    return series
end

"""
    modifying_integrator(rk::RungeKuttaMethod, order)
    modifying_integrator(A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                         order)

Compute the B-series of the [`modifying_integrator`](@ref) equation of the
Runge-Kutta method with Butcher coefficients `A, b, c` up to the prescribed
`order`.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be
    multiplied by a power of the time step divided by the `symmetry` of the
    rooted tree and multiplied by the corresponding elementary differential
    of the input vector field ``f``.
    See also [`evaluate`](@ref).
"""
function modifying_integrator(rk::RungeKuttaMethod, order)
    # B-series of the Runge-Kutta method
    series = bseries(rk, order)
    modifying_integrator(series)
end

function modifying_integrator(A::AbstractMatrix, b::AbstractVector,
                              c::AbstractVector, order)
    rk = RungeKuttaMethod(A, b, c)
    modifying_integrator(rk, order)
end

"""
    modifying_integrator(f, u, dt, series_integrator)

Compute the B-series of the [`modifying_integrator`](@ref) equation of the
time integration method with B-series `series_integrator` with respect to the
ordinary differential equation ``u'(t) = f(u(t))`` with vector field `f` and
dependent variables `u` for a time step size `dt`.

Here, `u` is assumed to be a vector of symbolic variables and `f` is assumed
to be a vector of expressions in these variables for plain B-series. For
B-series with colored trees, `f` must be a tuple of vectors of expressions in
the variables `u`. Currently, symbolic variables from

- [SymEngine.jl](https://github.com/symengine/SymEngine.jl),
- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), and
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

are supported.
"""
function modifying_integrator(f, u, dt, series_integrator)
    series = modifying_integrator(series_integrator)
    reduce_order_by = 1 # reduce the powers of `dt` by one since dt*f is given
    # by recursively solving `substitute`
    evaluate(f, u, dt, series, reduce_order_by)
end

"""
    modifying_integrator(f, u, dt, rk::RungeKuttaMethod, order)
    modifying_integrator(f, u, dt,
                         A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                         order)

Compute the B-series of the [`modifying_integrator`](@ref) equation of the
Runge-Kutta method with Butcher coefficients `A, b, c` up to the prescribed
`order` with respect to the ordinary differential equation ``u'(t) = f(u(t))``
with vector field `f` and dependent variables `u` for a time step size `dt`.

Here, `u` is assumed to be a vector of symbolic variables and `f` is assumed
to be a vector of expressions in these variables for plain B-series. For
B-series with colored trees, `f` must be a tuple of vectors of expressions in
the variables `u`. Currently, symbolic variables from

- [SymEngine.jl](https://github.com/symengine/SymEngine.jl),
- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), and
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

are supported.
"""
function modifying_integrator(f, u, dt, rk::RungeKuttaMethod, order)
    series_integrator = bseries(rk, order)
    modifying_integrator(f, u, dt, series_integrator)
end

function modifying_integrator(f, u, dt,
                              A::AbstractMatrix, b::AbstractVector,
                              c::AbstractVector, order)
    rk = RungeKuttaMethod(A, b, c)
    modifying_integrator(f, u, dt, rk, order)
end

"""
    elementary_differentials(f::AbstractVector, u, order)

Compute all elementary differentials of the vector field `f` with independent
variables `u` up to the given `order`. The return value can be indexed by
rooted trees to obtain the corresponding elementary differential.
"""
function elementary_differentials(f::AbstractVector, u, order)
    order >= 1 || throw(ArgumentError("The `order` must be at least one (got $order)"))
    differentials = OrderedDict{RootedTree{Int, Vector{Int}}, typeof(f)}()

    t = rootedtree(Int[])
    differentials[t] = one.(f)

    t = rootedtree([1])
    differentials[t] = f

    # Compute all necessary partial derivatives at first
    derivatives = Array{eltype(f)}[]
    push!(derivatives, f)
    for o in 1:(order - 1)
        d = similar(f, eltype(f), (length(f), ntuple(_ -> length(u), o)...))
        _compute_partial_derivatives!(d, f, u, derivatives[o])
        push!(derivatives, d)
    end

    # Compute all elementary differentials
    for o in 2:order
        for _t in RootedTreeIterator(o)
            t = copy(_t)
            differentials[t] = elementary_differential(f, t, differentials, derivatives)
        end
    end

    return differentials
end

function _compute_partial_derivatives!(d, f, u, lower_derivatives)
    for idx in CartesianIndices(d)
        idx_tuple = Tuple(idx)
        u_idx = Base.tail(idx_tuple)

        # All this B-series analysis only really makes sense for smooth functions.
        # Hence, we can use the symmetry of the partial derivatives to speed-up
        # the computations - Hermann Amandus Schwarz helps us again!
        issorted(u_idx) || continue

        # Next, we re-use already computed partial derivatives. Thus, we do not
        # compute the full set of derivatives as in
        #   f_idx = first(idx_tuple)
        #   partial_derivative = f[f_idx]
        #   for i in u_idx
        #     partial_derivative = compute_derivative(partial_derivative, u[i])
        #   end
        # but we re-use the already computed `lower_derivatives`.
        idx_known = Base.front(idx_tuple)
        idx_new = last(idx_tuple)
        partial_derivative = compute_derivative(lower_derivatives[idx_known...], u[idx_new])

        d[idx] = partial_derivative
    end
end

"""
    compute_derivative(expression, variable)

Internal function specialized on symbolic variables and expressions from

- [SymEngine.jl](https://github.com/symengine/SymEngine.jl),
- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), and
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

if these packages are loaded (via Requires.jl or weak dependencies on
Julia v1.9 and newer).
"""
function compute_derivative end

@static if !isdefined(Base, :get_extension)
    function __init__()
        @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" begin
            include("../ext/SymbolicsExt.jl")
        end

        @require SymEngine="123dc426-2d89-5057-bbad-38513e3affd8" begin
            include("../ext/SymEngineExt.jl")
        end

        @require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" begin
            include("../ext/SymPyExt.jl")
        end
    end
end

function elementary_differential(f, t, differentials, derivatives)
    subtr = RootedTrees.subtrees(t)
    result = similar(f)

    input_differentials = ntuple(n -> differentials[subtr[n]], length(subtr))
    input_derivative = derivatives[length(subtr) + 1]

    _compute_elementary_differential!(result, input_differentials, input_derivative)

    return result
end

@generated function _compute_elementary_differential!(result, input_differentials,
                                                      input_derivative::AbstractArray{T, N}) where {
                                                                                                    T,
                                                                                                    N
                                                                                                    }
    quote
        for i in eachindex(result)
            val = zero(eltype(result))

            Base.Cartesian.@nloops $(N - 1) j input_derivative begin
                # A naive version uses the full input of all possible combinations of
                # partial derivatives. If these were all computed in
                # `_compute_partial_derivatives!` above, we could just use
                #   term = Base.Cartesian.@ncall $(N-1) getindex input_derivative i j
                # and get rid of `sorted_j` here. However, all this B-series analysis
                # only really makes sense for smooth functions. Hence, we can use the
                # symmetry of the partial derivatives to speed-up the computations
                # - Hermann Amandus Schwarz helps us again!
                # TODO: Shall we use the non-allocating version `TupleTools.sort∘tuple`
                # instead of `sort!∘Base.vect`?
                sorted_j = Base.Cartesian.@ncall $(N - 1) sort!∘Base.vect j
                term = Base.Cartesian.@ncall $(N - 1) getindex input_derivative i d->sorted_j[d]
                Base.Cartesian.@nexprs $(N - 1) d->term *= input_differentials[d][j_d]
                val += term
            end

            result[i] = val
        end
    end
end

"""
    elementary_differentials(fs::NTuple{2, AbstractVector}, u, order)

Compute all elementary differentials of the sum of the two vector fields `f`
with independent variables `u` up to the given `order`. The return value
can be indexed by (bi-) colored rooted trees to obtain the corresponding
elementary differential.
"""
function elementary_differentials(f::NTuple{2, AbstractVector}, u, order)
    order >= 1 || throw(ArgumentError("The `order` must be at least one (got $order)"))
    N = 2 # length(fs)

    # Empty bicolored tree
    t = rootedtree(Int[], Bool[])
    differentials = OrderedDict{typeof(t), promote_type(map(typeof, f)...)}()
    differentials[t] = one.(f[1])

    # Bicolored trees with a single node
    t = rootedtree([1], Bool[0])
    differentials[t] = f[1]

    t = rootedtree([1], Bool[1])
    differentials[t] = f[2]

    # Compute all necessary partial derivatives at first
    derivatives = ntuple(n -> Array{eltype(f[n])}[], N)
    for n in 1:N
        push!(derivatives[n], f[n])
    end
    for o in 1:(order - 1)
        for n in 1:N
            d = similar(f[n], eltype(f[n]), (length(f[n]), ntuple(_ -> length(u), o)...))
            _compute_partial_derivatives!(d, f[n], u, derivatives[n][o])
            push!(derivatives[n], d)
        end
    end

    # Compute all elementary differentials
    for o in 2:order
        for _t in BicoloredRootedTreeIterator(o)
            t = copy(_t)
            differentials[t] = elementary_differential(f, t, differentials, derivatives)
        end
    end

    return differentials
end

function elementary_differential(f::NTuple{2, AbstractVector},
                                 t::ColoredRootedTree,
                                 differentials, derivatives)
    subtr = RootedTrees.subtrees(t)
    n = RootedTrees.color_to_index(root_color(t))
    result = similar(f[n])

    input_differentials = ntuple(i -> differentials[subtr[i]], length(subtr))
    input_derivative = derivatives[n][length(subtr) + 1]

    _compute_elementary_differential!(result, input_differentials, input_derivative)

    return result
end

include("latexify.jl")

# explicit precompilation on Julia v1.8 and newer
@static if VERSION >= v"1.8"
    include("precompile.jl")
end

"""
    renormalize!(series)

This function modifies a B-series by dividing each coefficient
by the `symmetry` of the corresponding tree.

!!! warning
    This breaks assumptions made on the representation of a B-series.
    The modified B-series should not be passed to any other function assuming
    the default normalization.
"""
function renormalize!(series)
    for t in keys(series)
        series[t] /= symmetry(t)
    end
    return series
end

"""
    energy_preserving_order(rk::RungeKuttaMethod, max_order; tol::Real=0)

This function checks up to which order a Runge-Kutta method `rk` is
energy-preserving for Hamiltonian problems.
It requires a `max_order` so that it does not run forever if the order up to
which the method is energy-preserving is too big or infinite.
Keyword Arguments:
- `tol::Real=1e-14`: The absolute tolerance for energy preservation.

See also [`is_energy_preserving`](@ref)
"""
function energy_preserving_order(rk::RungeKuttaMethod, max_order; tol::Real=0)
    p = 0
    not_energy_preserving = false
    while not_energy_preserving == false
        # Make sure not to pass the max_order
        if p > max_order
            return max_order
        end
        # Check energy preservation up to the given order
        if is_energy_preserving(rk, p + 1, tol=tol) == false
            not_energy_preserving = true
        end
        p = p + 1
    end
    return p - 1
end

"""
    is_energy_preserving(rk::RungeKuttaMethod, order; tol::Real=0)::Bool

This function checks whether the Runge-Kutta method `rk` is
energy-preserving for Hamiltonian systems up to a given `order`.
Keyword Arguments:
- `tol::Real=1e-14`: The absolute tolerance for energy preservation.
"""
function is_energy_preserving(rk::RungeKuttaMethod, order; tol::Real=0)
    series = bseries(rk, order)
    return is_energy_preserving(series, tol=tol)
end

"""
    is_energy_preserving(series_integrator; tol::Real=0)::Bool

This function checks whether the B-series `series_integrator` of a time
integration method is energy-preserving for Hamiltonian systems - up to the
[`order`](@ref) of `series_integrator`.
Keyword Arguments:
- `tol::Real=1e-14`: The absolute tolerance for energy preservation.

# References

This code is based on the Theorem 2 of
- Elena Celledoni, Robert I. McLachlan, Brynjulf Owren, and G. R. W. Quispel.
  "Energy-preserving integrators and the structure of B-series."
  Foundations of Computational Mathematics 10 (2010): 673-693.
  [DOI: 10.1007/s10208-010-9073-1](https://link.springer.com/article/10.1007/s10208-010-9073-1)
"""
function is_energy_preserving(series_integrator; tol::Real=0)
    # This method requires the B-series of a map as input
    let t = first(keys(series_integrator))
        @assert isempty(t)
        @assert isone(series_integrator[t])
    end

    # Theorem 2 of Celledoni et al. (2010) requires working with the modified
    # equation. The basic idea is to check whether the modified equation lies
    # in a certain subspace spanned by energy-preserving pairs of trees.
    # This implementation checks this condition for each order.
    flow_series = modified_equation(series_integrator)

    # It is easier to work with renormalized coefficients such that we do not
    # have to divide by the symmetry of the trees in the following.
    renormalize!(flow_series)

    # Save all the coefficients in an array: it is easier to use an array
    # of coefficients in '_is_energy_preserving'
    coefficients = collect(values(flow_series))
    # Save all the RootedTrees in another array
    thetrees = collect(keys(flow_series))
    # We need only the level sequences
    # Create an empty vector to store the converted trees into arrays
    trees = similar(thetrees, typeof(thetrees[1].level_sequence))
    # Convert the trees and store them in the `trees` vector
    for i in eachindex(trees, thetrees)
        trees[i] = thetrees[i].level_sequence
    end
    max_length = length(trees[end])

    # The `rank` function used in `_is_energy_preserving` doesn't work individually
    # for order less than 4 because the energy-preserving basis has only one element.
    # The following function checks the energy-preserving condition up to order 4
    # together.
    if _is_energy_preserving_low_order(trees, coefficients, tol) == false
        return false
    end

    # Next, we check each higher order separately. This reduces the problem
    # from one big problem to several smaller problems.
    for order in 5:max_length
        same_order_trees = empty(trees)
        same_order_coeffs = empty(coefficients)
        # We create and fill the arrays of trees and coefficients
        # corresponding to an `order`
        for j in eachindex(trees, coefficients)
            if length(trees[j]) == order
                push!(same_order_coeffs, coefficients[j])
                push!(same_order_trees, trees[j])
            end
        end
        if _is_energy_preserving(same_order_trees, same_order_coeffs, tol) == false
            return false
        end
    end

    return true
end

# This function checks whether the set of level_sequences `trees`
# and the set of 'coefficients' corresponding to a B-series is
# energy_preserving up to the 'uppermost_order', which we choose
# to be 4 for computational optimization
function _is_energy_preserving_low_order(trees, coefficients, tol)
    # Set the limit up to which this function will work
    uppermost_order = 4
    same_order_trees = empty(trees)
    same_order_coeffs = empty(coefficients)
    for index in eachindex(trees, coefficients)
        t = trees[index]
        if length(t) > uppermost_order
            break
        end
        push!(same_order_trees, t)
        push!(same_order_coeffs, coefficients[index])
    end
    return _is_energy_preserving(same_order_trees, same_order_coeffs, tol)
end

# This function is the application of Theorem 2 of the paper
# "Energy-Preserving Integrators and the Structure of B-series".
# It takes an array `trees` of `level_sequence`s and the array
# of coefficients.
# Checks whether `trees` and `ocefficients` satisfify the energy_preserving
# condition.
function _is_energy_preserving(trees, coefficients, tol)
    # TODO: `Float32` would also be nice to have. However, the default tolerance
    #       of the rank computation is based on `Float64`. Thus, it will usually
    #       not work with coefficients given only in 32 bit precision.
    # For very few coefficients, dense operations are more efficient than
    # sparse operations.
    if (length(coefficients) > 20 &&
        eltype(coefficients) <: Union{Float64,
              Rational{Int8}, Rational{Int16}, Rational{Int32}, Rational{Int64},
              Rational{Int128}})
        # These types support efficient computations in sparse matrices
        _is_energy_preserving_sparse(trees, coefficients, tol)
    else
        _is_energy_preserving_dense(trees, coefficients, tol)
    end
end

function _is_energy_preserving_dense(trees, coefficients, tol)
    # For every tree, obtain the adjoint and check if it exists
    length_coeff = length(trees)
    # Provided that a tree can have several adjoints, for every tree `t`
    # we need to check if there is a linear combination of the coefficients
    # of its adjoints such that this linear combination is equal to the
    # coefficient of `t`.
    # For this purpose, we create a energy_preserving_basis to store all the
    # information
    energy_preserving_basis = empty(trees)
    for (t_index, t) in enumerate(trees)
        # We do not need to consider the empty tree
        isempty(t) && continue

        # Check whether the level sequence corresponds to a bushy tree.
        # If so, we can exit early since the coefficients of bushy trees
        # must be zero in the modified equation of energy-preserving methods.
        if length(t) > 1 && is_bushy(t)
            if (typeof(coefficients[t_index]) == Float64 || typeof(coefficients[t_index]) == Float32 || typeof(coefficients[t_index]) == Float16) && tol == 0
                tol = 1e-14
                if !isapprox(coefficients[t_index], 0, atol = tol)
                    return false
                end
            elseif tol != 0
                if !isapprox(coefficients[t_index], 0, atol = tol)
                    return false
                end
            else
                if !iszero(coefficients[t_index])
                    return false
                end
            end
        else
            # If the tree is not a bush, then generate all the equivalent trees
            equiv_set = equivalent_trees(t)
            for onetree in equiv_set
                # j-th canonical vector
                ej = zeros(Int64, length_coeff)
                ej[t_index] = 1
                m = num_branches(onetree)
                energy_preserving_pair = rightmost_energy_preserving_tree(onetree)
                ek = zeros(Int64, length_coeff)
                # Check if the rightmost energy-preserving tree is in the
                # set of trees
                idx = findfirst(isequal(energy_preserving_pair), trees)
                if idx !== nothing # i.e., energy_preserving_pair in trees
                    # Generate a k-th canonical vector
                    ek[idx] = 1 * (-1)^m
                    # Add ej + ek and push it into energy_preserving_basis
                    ej = ej + ek
                    push!(energy_preserving_basis, ej)
                end
            end
        end
    end
    # We remove the empty columns (those full of zeros)
    filter!(!iszero, energy_preserving_basis)

    # We remove repeated columns
    unique!(energy_preserving_basis)

    # The components of `energy_preserving_basis` are the columns of the matrix `M`.
    M = reduce(hcat, energy_preserving_basis)
    rank_M = rank(M)

    # We also create an extended matrix for which we append the vector of
    # coefficients
    Mv = [M coefficients]
    rank_Mv = rank(Mv)

    # If the rank of M is equal to the rank of the extended Mv,
    # then the system is energy-preserving
    return rank_M == rank_Mv
end

# This method is specialized for coefficients that can be used efficiently in
# sparse arrays.
function _is_energy_preserving_sparse(trees, coefficients, tol)
    # For every tree, obtain the adjoint and check if it exists.
    # Provided that a tree can have several adjoints, for every tree `t`
    # we need to check if there is a linear combination of the coefficients
    # of its adjoints such that this linear combination is equal to the
    # coefficient of `t`.
    # For this purpose, we create a energy-preserving basis as sparse matrix
    # to store all the information. We assemble this matrix in a sparse format.
    rows = Vector{Int}()
    cols = Vector{Int}()
    # The SparseArrays package supports only floating point coefficients.
    vals = Vector{Float64}()

    col = 1
    for (t_index, t) in enumerate(trees)
        # We do not need to consider the empty tree
        isempty(t) && continue

        # Check whether the level sequence corresponds to a bushy tree.
        # If so, we can exit early since the coefficients of bushy trees
        # must be zero in the modified equation of energy-preserving methods.
        if length(t) > 1 && is_bushy(t)
            if (typeof(coefficients[t_index]) == Float64 || typeof(coefficients[t_index]) == Float32 || typeof(coefficients[t_index]) == Float16) && tol == 0
                tol = 1e-14
                if !isapprox(coefficients[t_index], 0, atol = tol)
                    return false
                end
            elseif tol != 0
                if !isapprox(coefficients[t_index], 0, atol = tol)
                    return false
                end
            else
                if !iszero(coefficients[t_index])
                    return false
                end
            end
        else
            # If the tree is not a bush, then generate all the equivalent trees
            equiv_set = equivalent_trees(t)
            for onetree in equiv_set
                push!(rows, t_index)
                push!(cols, col)
                push!(vals, 1)

                m = num_branches(onetree)
                energy_preserving_pair = rightmost_energy_preserving_tree(onetree)

                # Check if the rightmost energy-preserving tree is in the
                # set of trees
                idx = findfirst(isequal(energy_preserving_pair), trees)
                if idx !== nothing # i.e., energy_preserving_pair in trees
                    push!(rows, idx)
                    push!(cols, col)
                    push!(vals, (-1)^m)
                    col = col + 1
                end
            end
        end
    end
    # TODO: We could try to improve the performance by
    # - removing empty columns (those full of zeros)
    # - removing repeated columns
    # - working with `svd` directly instead of computing the `rank` twice

    # The components of `energy_preserving_basis` are the columns of the matrix `M`.
    M = sparse(rows, cols, vals, length(coefficients), maximum(cols))
    rank_M = rank(M)

    # We also create an extended matrix for which we append the vector of
    # coefficients
    col = maximum(cols) + 1
    for i in eachindex(coefficients)
        coeff = coefficients[i]
        iszero(coeff) && continue
        push!(rows, i)
        push!(cols, col)
        push!(vals, coeff)
    end
    Mv = sparse(rows, cols, vals, length(coefficients), col)
    rank_Mv = rank(Mv)

    # If the rank of M is equal to the rank of the extended Mv,
    # then the system is energy-preserving
    return rank_M == rank_Mv
end

# This function checks if the level_sequence is a bush.
function is_bushy(tree)
    return maximum(tree) <= tree[1] + 1
end

function num_branches(a)
    k = a[end]
    return k - 1
end

# This function returns the rightmost energy-preserving tree (`level_sequence`)
# for a given tree (the input also in level-sequence form).
function rightmost_energy_preserving_tree(a::Vector{Int})
    # We want to generate all the branches with respect to the rightmost trunk
    branches = get_branches(a)
    # We obtain the number of branches the right-most trunk has.
    # Theorem 2 in the article requires to know if m is odd or even
    m = num_branches(a)
    # We create an array from 1 to m plus another node for the final branch
    # of the rightmost trunk
    energy_preserving_partner = collect(1:(m + 1))
    # Then, we insert every level_sequence from branches
    for j in 1:m
        last_j_occurrence::Int = findlast(x -> x == j, energy_preserving_partner)
        last_jplus1_occurrence::Int = findlast(x -> x == j + 1, energy_preserving_partner)
        energy_preserving_partner::typeof(a) = vcat(energy_preserving_partner[1:last_j_occurrence],
                                                    branches[j],
                                                    energy_preserving_partner[last_jplus1_occurrence:end])
    end
    energy_preserving_partner = canonicalarray(energy_preserving_partner)
    return energy_preserving_partner
end

# This function returns the forest of `branches` after removing the ritghtmost-trunk
# from the tree. This array also swaps the order of the branches and modifies
# the numbers in its level_sequence so that the new branches are ready for being
# reassembled.
function get_branches(a)
    m = num_branches(a)
    branches = Array{Array}(undef, m)
    # We need to save the final number in the level_sequence because this is
    # the final branch of the trunk
    k = a[end]
    # We need to look for the last `last_j_occurrence` of every integer in [1, k-1]
    for j in 1:(k - 1)
        last_j_occurrence = findlast(x -> x == j, a)
        last_jplus1_occurrence = findlast(x -> x == j + 1, a)
        # Consider the empty branches
        if isnothing(last_j_occurrence) || isnothing(last_jplus1_occurrence)
            branches[j] = []
        else
            branches[j] = a[(last_j_occurrence + 1):(last_jplus1_occurrence - 1)]
        end
    end
    # Create another dict for the modified indexes
    modified_branches = Array{Array}(undef, m)
    mid_branch = (m + 1) / 2
    # We check if the number of branches is odd:
    # in that case, the middle one remains the same
    for j in 1:m
        if m % 2 == 1 && j == mid_branch
            modified_branches[j] = branches[j]
            # We use the formula
            # n+m-2j+1 for every number in the level_sequence
        else
            modified_branches[j] = [n + m - 2j + 1 for n in branches[j]]
        end
    end
    return reverse(modified_branches)
end

# This function rearranges the level sequence in
# the same way as the BSeries.jl output does
function canonicalarray(tree)
    t = rootedtree(tree)
    trees = t.level_sequence
    return trees
end

# This function generates all the equivalent trees for a given `level_sequence`
function equivalent_trees(tree)
    equivalent_trees_set = Vector{typeof(tree)}()
    l = length(tree)
    for i in 1:(l - 1)
        # For every node check if it is a leaf
        if tree[i + 1] <= tree[i]
            # Make a copy of `tree` to work with
            graft = copy(tree)
            leaf_value = graft[i]
            # Initialize the length of the subtree we are currently working with
            length_subtree = l
            # Initialize the initial component from which we will start looping
            # in order to cut and paste every subtree at the right of the array
            initial_component = 1
            left_index = 0
            right_index = 0
            # Sending subtrees to the right
            leaf_current_position = i
            for j in 2:leaf_value
                # Make sure that there are more than one 'j'
                # in the current level sequence
                j_counter = 0
                for node in initial_component:l
                    if graft[node] == j
                        j_counter += 1
                    end
                end
                # If there is only one 'j', then it makes no sense
                # to continue the loop: continue with the next 'j'
                if j_counter == 1
                    continue
                end
                # Find the nearest 'j' on the left side, from 'initial_component'
                # until 'leaf_current_position'
                for k in leaf_current_position:-1:initial_component
                    if graft[k] == j
                        left_index = k
                        break
                    end
                end
                # Find the nearest 'j' on the right side
                right_flag = false
                for k in leaf_current_position:l
                    if graft[k] == j
                        right_flag = true
                        right_index = k
                        break
                    end
                end
                # Special case
                if right_index == 0
                    continue
                end
                # Consider special case: because we are defining (below) the right limit
                # for cutting the tree to be (`right_index` - 1), if the leaf is near to
                # right border we must define the right_index to be l+1.
                if right_flag == false
                    right_index = l + 1
                end
                # Get the subtree to be transplanted
                if left_index == right_index
                    plantout_subtree = [j]
                    # Remove these components from the original tree
                    splice!(graft, left_index)
                else
                    plantout_subtree = graft[left_index:(right_index - 1)]
                    # Remove these components from the original tree
                    splice!(graft, left_index:(right_index - 1))
                end
                length_subtree = length(plantout_subtree)
                # Re-define the initial_component from which the subtree
                # we are working with starts
                initial_component = l - length_subtree + 1
                # Add the components to the right
                graft = vcat(graft, plantout_subtree)
                # Refresh the current position of the leaf
                leaf_current_position = l + leaf_current_position - length_subtree -
                                        left_index + 1
                if leaf_current_position == l
                    break
                end
            end
            # Obtain the energy_preserving_partner of 'graft'
            push!(equivalent_trees_set, graft)
        end
    end
    # Include the original tree
    push!(equivalent_trees_set, tree)
    # Eliminate repeated
    unique!(equivalent_trees_set)

    return equivalent_trees_set
end

end # module
