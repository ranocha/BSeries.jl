module BSeries

using Reexport: @reexport

@reexport using RootedTrees
using RootedTrees: RootedTree

@reexport using OrderedCollections: OrderedDict

using Requires: @require


export TruncatedBSeries, ExactSolution

export bseries, substitute, compose

export modified_equation, modifying_integrator

export elementary_differentials



# types used for traits
abstract type AbstractEvaluation end
struct LazyEvaluation     <: AbstractEvaluation end # lazy evaluation up to arbitrary order
struct EagerEvaluation    <: AbstractEvaluation end # eager evaluation up to a fixed order
struct MemoizedEvaluation <: AbstractEvaluation end # lazy evaluation up to arbitrary order memoizing the results

# internal interface
# TODO: `bseries(series::AbstractBSeries) = series` or a `copy`?
# TODO: `bseries(series, order)`` returns a `copy`? up to `order`, if available
# TODO: `truncate!(series, order)`


"""
    TruncatedBSeries

An [`AbstractBSeries`](@ref) that can describe B-series of both numerical
integration methods (where the coefficient of the empty tree is unity)
and right-hand sides of an ordinary differential equation and perturbations
thereof (where the coefficient of the empty tree is zero) up to a prescribed
[`order`](@ref).

Generally, it should be constructed via [`bseries`](@ref) or one of the
methods returning a B-series.
"""
struct TruncatedBSeries{T<:RootedTree, V} <: AbstractDict{T, V}
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

@inline Base.getindex(series::TruncatedBSeries, t::RootedTree) = getindex(series.coef, t)
@inline Base.setindex!(series::TruncatedBSeries, val, t::RootedTree) = setindex!(series.coef, val, t)

@inline Base.get(series::TruncatedBSeries, t::RootedTree, default) = get(series.coef, t, default)
@inline Base.get(f::F, series::TruncatedBSeries, t::RootedTree) where {F<:Function} = get(f, series.coef, t)

@inline Base.getkey(series::TruncatedBSeries, t::RootedTree, default) = getkey(series.coef, t, default)

@inline Base.delete!(series::TruncatedBSeries, t::RootedTree) = (delete!(series.coef, t); series)

@inline Base.pop!(series::TruncatedBSeries, t::RootedTree) = pop!(series.coef, t)
@inline Base.pop!(series::TruncatedBSeries, t::RootedTree, default) = pop!(series.coef, t, default)

@inline Base.sizehint!(series::TruncatedBSeries, n) = sizehint!(series.coef, n)


# internal interface of B-series
@inline evaluation_type(::TruncatedBSeries) = EagerEvaluation()

# TODO: TruncatedBSeries
# This assumes that the `coef` are always stored with increasing `order` of the
# rooted trees and that the underlying data format in OrderedCollections.jl does
# not change. A general fallback would be
#   RootedTrees.order(series::TruncatedBSeries) = maximum(order, keys(series.coef))
# but that is O(n) instead of O(1).
RootedTrees.order(series::TruncatedBSeries) = order(series.coef.keys[end])



"""
    ExactSolution{V}()

Lazy representation of the B-series of the exact solution of an ordinary
differential equation using coefficients of type at least as representative as
`V`.
"""
struct ExactSolution{V} end

Base.getindex(::ExactSolution{V}, t::RootedTree) where {V} = inv(convert(V, γ(t)))

# internal interface of B-series
@inline evaluation_type(::ExactSolution) = LazyEvaluation()

# general interface methods of iterators for `ExactSolution`
Base.IteratorSize(::Type{<:ExactSolution}) = Base.SizeUnknown()
Base.eltype(::Type{ExactSolution{V}}) where {V} = V

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



"""
    bseries(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)

Compute the B-series of the Runge-Kutta method with Butcher coefficients
`A, b, c` up to a prescribed `order`
"""
function bseries(A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                 order)
  V = promote_type(eltype(A), eltype(b), eltype(c))
  series = TruncatedBSeries{RootedTree{Int, Vector{Int}}, V}()

  for o in 1:order
    for t in RootedTreeIterator(o)
      series[copy(t)] = elementary_weight(t, A, b, c)
    end
  end

  return series
end

# TODO: bseries(A::AbstractMatrix, b::AbstractVector, c::AbstractVector)
# should create a lazy version, optionally a meoized one



"""
    substitute(b, a, t::RootedTree)

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
function substitute(b, a, t::RootedTree)
  result = zero(first(values(a)) * first(values(b)))

  for (forest, skeleton) in PartitionIterator(t)
    update = a[skeleton]
    for tree in forest
      update *= b[tree]
    end
    result += update
  end

  return result
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
    for tree in forest
      update *= b[tree]
    end
    result += update
  end

  return result
end



"""
    modified_equation(series_integrator)

Compute the B-series of the modified equation of the time integration method
with B-series `series_integrator`.

Given an ordinary differential equation (ODE) ``u'(t) = f(u(t))`` and a
Runge-Kutta method, the idea is to interpret the numerical solution with
given time step size as exact solution of a modified ODE ``u'(t) = fₕ(u(t))``.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be
    multiplied by the corresponding elementary differential of the input
    vector field ``f``.

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
modified_equation(series_integrator) = modified_equation(series_integrator,
                                                         evaluation_type(series_integrator))

function modified_equation(series_integrator, ::EagerEvaluation)
  # B-series of the exact solution
  V = valtype(series_integrator)

  # We could just use
  #   series_ex = ExactSolution{V}()
  # However, we need to access elements of `series_ex` more than once in the
  # subsitution below. Thus, it's cheaper to compute every entry only once and
  # re-use it later.
  exact = ExactSolution{V}()
  series_ex = empty(series_integrator)
  for t in keys(series_integrator)
    series_ex[t] = exact[t]
  end

  # Prepare B-series of the modified equation
  series = empty(series_integrator)
  for t in keys(series_integrator)
    series[t] = zero(V)
  end

  t = first(keys(series_integrator))
  series[t] = series_integrator[t]

  # Recursively solve `substitute(series, series_ex, t) == series_integrator[t]`.
  # This works because
  #   substitute(series, series_ex, t) = series[t] + lower order terms
  # Since the `keys` are ordered, we don't need to use nested loops of the form
  #   for o in 2:order
  #     for _t in RootedTreeIterator(o)
  #       t = copy(_t)
  # which are slightly less efficient due to additional computations and allocations.
  for t in keys(series_integrator)
    series[t] += series_integrator[t] - substitute(series, series_ex, t)
  end

  return series
end

"""
    modified_equation(A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                      order)

Compute the B-series of the [`modified_equation`](@ref) of the Runge-Kutta
method with Butcher coefficients `A, b, c` up to the prescribed `order`.
"""
function modified_equation(A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                           order)
  # B-series of the Runge-Kutta method
  series = bseries(A, b, c, order)
  modified_equation(series)
end


"""
    modified_equation(f, u, dt, series_integrator)

Compute the [`modified_equation`](@ref) of the time integration method with
B-series `series_integrator` specialized to the ordinary differential equation
``u'(t) = f(u(t))`` with vector field `f` and dependent variables `u` for a
time step size `dt`.

Here, `u` is assumed to be a vector of symbolic variables and `f` is assumed
to be a vector of expressions in these variables. Currently, symbolic variables
from

- [SymEngine.jl](https://github.com/symengine/SymEngine.jl),
- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), and
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

are supported.

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
modified_equation(f, u, dt, series_integrator)

function modified_equation(f, u, dt,
                           A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                           order)
  series = modified_equation(A, b, c, order)
  differentials = elementary_differentials(f, u, order)
  result = zero(f)
  for t in keys(series)
    result += dt^(RootedTrees.order(t) - 1) / symmetry(t) * series[t] * differentials[t]
  end
  result
end


"""
    modified_equation(f, u, dt,
                      A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)

Compute the B-series of the [`modified_equation`](@ref) of the Runge-Kutta
method with Butcher coefficients `A, b, c` up to the prescribed `order` with
respect to the ordinary differential equation ``u'(t) = f(u(t))`` with vector
field `f` and dependent variables `u` for a time step size `dt`.

Here, `u` is assumed to be a vector of symbolic variables and `f` is assumed
to be a vector of expressions in these variables. Currently, symbolic variables
from

- [SymEngine.jl](https://github.com/symengine/SymEngine.jl),
- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), and
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

are supported.

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function modified_equation(f, u, dt,
                           A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                           order)
  series = modified_equation(A, b, c, order)
  differentials = elementary_differentials(f, u, order)
  result = zero(f)
  for t in keys(series)
    result += dt^(RootedTrees.order(t) - 1) / symmetry(t) * series[t] * differentials[t]
  end
  result
end


"""
    modifying_integrator(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)

Compute the B-series of a "modifying integrator" equation of the Runge-Kutta
method with Butcher coefficients `A, b, c` up to the prescribed `order`.

Given an ordinary differential equation (ODE) ``u'(t) = f(u(t))`` and a
Runge-Kutta method, the idea is to find a modified ODE ``u'(t) = fₕ(u(t))``
such that the numerical solution with given time step size is the exact solution
of the original ODE.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be multiplied
    by the corresponding elementary differential of the input vector field ``f``.

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function modifying_integrator(A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                              order)
  # B-series of the Runge-Kutta method
  series_rk = bseries(A, b, c, order)

  # B-series of the exact solution
  V = valtype(series_rk)
  series_ex = ExactSolution{V}()

  # Prepare B-series of the modifying integrator equation
  series = empty(series_rk)
  for t in keys(series_rk)
    series[t] = zero(V)
  end

  t = rootedtree([1])
  series[t] = series_rk[t]

  # Recursively solve `substitute(series, series_rk, t) == series_ex[t]`.
  # This works because
  #   substitute(series, series_rk, t) = series[t] + lower order terms
  # Since the `keys` are ordered, we don't need to use nested loops of the form
  #   for o in 2:order
  #     for _t in RootedTreeIterator(o)
  #       t = copy(_t)
  # which are slightly less efficient due to additional computations and allocations.
  for t in keys(series)
    series[t] += series_ex[t] - substitute(series, series_rk, t)
  end

  return series
end


"""
    modifying_integrator(f, u, dt,
                         A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)

Compute the B-series of a "modifying integrator" equation of the Runge-Kutta
method with Butcher coefficients `A, b, c` up to the prescribed `order` with
respect to the ordinary differential equation ``u'(t) = f(u(t))`` with vector
field `f` and dependent variables `u` for a time step size `dt`.

Here, `u` is assumed to be a vector of symbolic variables and `f` is assumed
to be a vector of expressions in these variables. Currently, symbolic variables
from

- [SymEngine.jl](https://github.com/symengine/SymEngine.jl),
- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), and
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)

are supported.

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series.
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function modifying_integrator(f, u, dt,
                              A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                              order)
  series = modifying_integrator(A, b, c, order)
  differentials = elementary_differentials(f, u, order)
  result = zero(f)
  for t in keys(series)
    result += dt^(RootedTrees.order(t) - 1) / symmetry(t) * series[t] * differentials[t]
  end
  result
end



"""
    elementary_differentials(f, u, order)

Compute all elementary differentials of the vector field `f` with independent
variables `u` up to the given `order`. The return value can be indexed by
rooted trees to obtain the corresponding elementary differential.
"""
function elementary_differentials(f, u, order)
  order >= 1 || throw(ArgumentError("The `order` must be at least one (got $order)"))
  differentials = OrderedDict{RootedTree{Int, Vector{Int}}, typeof(f)}()

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
  end

  @require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" begin
    using .SymPy: SymPy

    function compute_derivative(expression::SymPy.Sym, variable::SymPy.Sym)
      SymPy.diff(expression, variable)
    end
  end

  @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" begin
    using .Symbolics: Symbolics

    function compute_derivative(expression::Symbolics.Num, variable::Symbolics.Num)
      Symbolics.expand_derivatives(Symbolics.Differential(variable)(expression))
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


end # module
