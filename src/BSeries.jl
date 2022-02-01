module BSeries

@doc read(joinpath(dirname(@__DIR__), "README.md"), String) BSeries


using Reexport: @reexport

@reexport using RootedTrees
using RootedTrees: AbstractRootedTree

@reexport using OrderedCollections: OrderedDict

using Requires: @require

using Latexify: Latexify, LaTeXString

@reexport using Polynomials: Polynomials, Polynomial


export TruncatedBSeries, ExactSolution

export bseries, substitute, compose, evaluate

export modified_equation, modifying_integrator

export elementary_differentials

export MultirateInfinitesimalSplitMethod



# Types used for traits
# These traits may decide between different algorithms based on the
# corresponding complexity etc.
abstract type AbstractEvaluation end
struct LazyEvaluation     <: AbstractEvaluation end # lazy evaluation up to arbitrary order
struct EagerEvaluation    <: AbstractEvaluation end # eager evaluation up to a fixed order
struct MemoizedEvaluation <: AbstractEvaluation end # lazy evaluation memoizing the results

evaluation_type(::Any) = LazyEvaluation() # assume lazy evaluation by default
evaluation_type(::AbstractDict) = EagerEvaluation() # dictionaries store results eagerly


# internal interface
# TODO: `bseries(series::AbstractBSeries) = series` or a `copy`?
# TODO: `bseries(series, order)`` returns a `copy`? up to `order`, if available
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
struct TruncatedBSeries{T<:AbstractRootedTree, V} <: AbstractDict{T, V}
  coef::OrderedDict{T, V}
end

TruncatedBSeries{T, V}() where {T, V} = TruncatedBSeries{T, V}(OrderedDict{T, V}())


# general interface methods of `AbstractDict` for `TruncatedBSeries`
@inline Base.iterate(series::TruncatedBSeries) = iterate(series.coef)
@inline Base.iterate(series::TruncatedBSeries, state) = iterate(series.coef, state)

@inline Base.empty(series::TruncatedBSeries) = TruncatedBSeries(empty(series.coef))
@inline Base.empty(series::TruncatedBSeries, ::Type{T}, ::Type{V}) where {T, V} = TruncatedBSeries(empty(series.coef, T, V))
@inline Base.empty!(series::TruncatedBSeries) = (empty!(series.coef); series)

@inline Base.length(series::TruncatedBSeries) = length(series.coef)

@inline Base.getindex(series::TruncatedBSeries, t::AbstractRootedTree) = getindex(series.coef, t)
@inline Base.setindex!(series::TruncatedBSeries, val, t::AbstractRootedTree) = setindex!(series.coef, val, t)

@inline Base.get(series::TruncatedBSeries, t::AbstractRootedTree, default) = get(series.coef, t, default)
@inline Base.get(f::F, series::TruncatedBSeries, t::AbstractRootedTree) where {F<:Function} = get(f, series.coef, t)

@inline Base.getkey(series::TruncatedBSeries, t::AbstractRootedTree, default) = getkey(series.coef, t, default)

@inline Base.delete!(series::TruncatedBSeries, t::AbstractRootedTree) = (delete!(series.coef, t); series)

@inline Base.pop!(series::TruncatedBSeries, t::AbstractRootedTree) = pop!(series.coef, t)
@inline Base.pop!(series::TruncatedBSeries, t::AbstractRootedTree, default) = pop!(series.coef, t, default)

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
"""
RootedTrees.order(series::TruncatedBSeries) = order(series.coef.keys[end])



"""
    ExactSolution{V}()

Lazy representation of the B-series of the exact solution of an ordinary
differential equation using coefficients of type at least as representative as
`V`.
"""
struct ExactSolution{V} end

Base.getindex(::ExactSolution{V}, t::AbstractRootedTree) where {V} = convert(V, 1//1) / γ(t)

# general interface methods of iterators for `ExactSolution`
Base.IteratorSize(::Type{<:ExactSolution}) = Base.SizeUnknown()
Base.eltype(::Type{ExactSolution{V}}) where {V} = V
Base.valtype(exact::ExactSolution) = valtype(typeof(exact))
Base.valtype(::Type{ExactSolution{V}}) where {V} = V

function Base.iterate(exact::ExactSolution)
  iterator = RootedTreeIterator(1)
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


# TODO: Documentation, Base.show, export etc.
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
struct MultirateInfinitesimalSplitMethod{T, PolyMatT<:AbstractMatrix{<:Polynomial{T}}, MatT<:AbstractMatrix{T}, VecT<:AbstractVector{T}} <: RootedTrees.AbstractTimeIntegrationMethod
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

  prototype_scalar = TruncatedBSeries{BicoloredRootedTree{Int, Vector{Int}, Vector{Bool}}, V}()
  prototype_scalar[rootedtree(Int[], Bool[])] = one(V)

  poly_one = Polynomial{V, :x}([one(V)])
  poly_zero = zero(poly_one)
  prototype_polynomial = TruncatedBSeries{BicoloredRootedTree{Int, Vector{Int}, Vector{Bool}}, typeof(poly_one)}()
  prototype_polynomial[rootedtree(Int[], Bool[])] = poly_one

  A = mis.A
  D = mis.D
  G = mis.G
  ns = size(A, 2)
  Z = Vector{typeof(prototype_polynomial)}(undef, ns + 1)
  η = Vector{typeof(prototype_scalar)}(undef, ns + 1)
  for i in 1:ns+1
    Z[i] = copy(prototype_polynomial)
    η[i] = copy(prototype_scalar)
  end

  for o in 1:order
    for t_ in BicoloredRootedTreeIterator(o)
      t = copy(t_)
      for i in 1:ns+1
        phi = poly_zero
        for j in 1:i-1
          r = poly_one
          for subtree in RootedTrees.SubtreeIterator(t)
            if subtree.color_sequence[1]
              r = r * Z[i][subtree]
            else
              r = r * η[j][subtree]
            end
          end

          v = [one(V)]
          for k in 0:length(A[i, j])-1
             phi = phi + Polynomial{V, :x}(v) * A[i, j][k] * r
             v = [0; v]
          end
        end

        for j in 1:i-1
          phi = phi + G[i,j] * (η[j][t] - η[1][t])
        end

        phi = Polynomials.integrate(phi)
        phi[0] = 0
        for j in 1:i-1
          phi[0] = phi[0] + D[i,j] * (η[j][t] - η[1][t])
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

Compute the coefficient correspoding to the tree `t` of the B-series that is
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

Compute the coefficient correspoding to the tree `t` of the B-series that is
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
`compose(b, a) = b ⋅ a`.

If `normalize_stepsize = true`, the coefficients of the returned B-series will
are divied by `2^order(t)` for each rooted tree `t`. This normalizes the step
size so that the resulting numerical integrator B-series uses the same step size
as the input series (instead of a doubled step size).

# References

Section 3.1 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function compose(b, a; normalize_stepsize=false)
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
where the B-series coeffcients are those of ``h fₕ``, i.e., they contain an
additional power of `dt`. In this case, the B-series of the vector field can
be obtained using `reduce_order_by = 1`.

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function evaluate(f, u, dt, series, reduce_order_by=0)
  evaluate(f, u, dt, series, evaluation_type(series), reduce_order_by)
end

function evaluate(f, u, dt, series, ::EagerEvaluation, reduce_order_by)
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
  modified_equation(series_integrator, evaluation_type(series_integrator))
end

function modified_equation(series_integrator, ::EagerEvaluation)
  V = valtype(series_integrator)

  # B-series of the exact solution
  # We could just use the lazy version
  #   series_ex = ExactSolution{V}()
  # However, we need to access elements of `series_ex` more than once in the
  # subsitution below. Thus, it's cheaper to compute every entry only once and
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
  modifying_integrator(series_integrator, evaluation_type(series_integrator))
end

function modifying_integrator(series_integrator, ::EagerEvaluation)
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
  for o in 1:(order-1)
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

if these packages are loaded (via Requires.jl).
"""
function compute_derivative end

function __init__()
  @require SymEngine="123dc426-2d89-5057-bbad-38513e3affd8" begin
    using .SymEngine: SymEngine

    function compute_derivative(expression::SymEngine.Basic, variable::SymEngine.Basic)
      SymEngine.diff(expression, variable)
    end

    latexify_default_dt(::Type{SymEngine.Basic}) = SymEngine.symbols("h")
  end

  @require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" begin
    using .SymPy: SymPy

    function compute_derivative(expression::SymPy.Sym, variable::SymPy.Sym)
      SymPy.diff(expression, variable)
    end

    latexify_default_dt(::Type{SymPy.Sym}) = SymPy.symbols("h", real=true)
  end

  @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" begin
    using .Symbolics: Symbolics

    function compute_derivative(expression::Symbolics.Num, variable::Symbolics.Num)
      Symbolics.expand_derivatives(Symbolics.Differential(variable)(expression))
    end

    # nested macro calls do not work so we need to define this part in another file
    include("latexify_symbolics.jl")
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

@generated function _compute_elementary_differential!(result, input_differentials, input_derivative::AbstractArray{T,N}) where {T,N}
  quote
    for i in eachindex(result)
      val = zero(eltype(result))

      Base.Cartesian.@nloops $(N-1) j input_derivative begin
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
        sorted_j = Base.Cartesian.@ncall $(N-1) sort!∘Base.vect j
        term = Base.Cartesian.@ncall $(N-1) getindex input_derivative i d -> sorted_j[d]
        Base.Cartesian.@nexprs $(N-1) d -> term *= input_differentials[d][j_d]
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
  for o in 1:(order-1)
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


end # module
