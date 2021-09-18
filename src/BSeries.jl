module BSeries

using Reexport: @reexport

@reexport using RootedTrees
using RootedTrees: RootedTree

@reexport using OrderedCollections: OrderedDict


export bseries, substitute

export modified_equation, modifying_integrator


"""
    substitute(b, a, t::RootedTree)

Compute the coefficient correspoding to the tree `t` of the B-series that is
formed by substituting the B-series `b` into the B-series `a`.
"""
function substitute(b, a, t::RootedTree)
  forests, skeletons = all_partitions(t)
  result = zero(first(values(a)) * first(values(b)))

  for (forest, skeleton) in zip(forests, skeletons)
    result += reduce(*, b[tree] for tree in forest) * a[skeleton]
  end

  return result
end



"""
    bseries(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)

Compute the B-series of the  Runge-Kutta method with Butcher coefficients
`A, b, c` up to a prescribed `order`
"""
function bseries(A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                 order)
  T = promote_type(eltype(A), eltype(b), eltype(c))
  series = OrderedDict{RootedTree{Int, Vector{Int}}, T}()

  for o in 1:order
    for t in RootedTreeIterator(o)
      series[copy(t)] = elementary_weight(t, A, b, c)
    end
  end

  return series
end


"""
    ExactSolution{T}

Lazy representation of the B-series of the exact solution of an ordinary
differential equation using coefficients of type at least as representative as
`T`.
"""
struct ExactSolution{T} end

Base.getindex(::ExactSolution{T}, t::RootedTree) where {T} = one(T) / γ(t)

Base.IteratorSize(::Type{<:ExactSolution}) = Base.SizeUnknown()
Base.eltype(::Type{ExactSolution{T}}) where {T} = promote_type(T, Int)

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
    modified_equation(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)

Compute the B-series of the modified equation of the Runge-Kutta method with
Butcher coefficients `A, b, c` up to the prescribed `order`.

Given an ordinary differential equation (ODE) ``u'(t) = f(u(t))`` and a
Runge-Kutta method, the idea is to interpret the numerical solution with
given time step size as exact solution of a modified ODE ``u'(t) = fₕ(u(t))``.

!!! note "Normalization by elementary differentials"
    The coefficients of the B-series returned by this method need to be multiplied
    by the corresponding elementary differential of the input vector field ``f``.

See Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function modified_equation(A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                           order)
  # B-series of the Runge-Kutta method
  series_rk = bseries(A, b, c, order)

  # B-series of the exact solution
  T = eltype(values(series_rk))
  series_ex = ExactSolution{T}()

  # Prepare B-series of the modified equation
  series = empty(series_rk)
  for t in keys(series_rk)
    series[t] = zero(T)
  end

  t = rootedtree([1])
  series[t] = series_rk[t]

  # Recursively solve `substitute(series, series_ex, t) == series_rk[t]`.
  # This works because
  #   substitute(series, series_ex, t) = series[t] + lower order terms
  for o in 2:order
    for _t in RootedTreeIterator(o)
      t = copy(_t)
      series[t] += series_rk[t] - substitute(series, series_ex, t)
    end
  end

  return series
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

See Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function modifying_integrator(A::AbstractMatrix, b::AbstractVector, c::AbstractVector,
                              order)
  # B-series of the Runge-Kutta method
  series_rk = bseries(A, b, c, order)

  # B-series of the exact solution
  T = eltype(values(series_rk))
  series_ex = ExactSolution{T}()

  # Prepare B-series of the modifying integrator equation
  series = empty(series_rk)
  for t in keys(series_rk)
    series[t] = zero(T)
  end

  t = rootedtree([1])
  series[t] = series_rk[t]

  # Recursively solve `substitute(series, series_rk, t) == series_ex[t]`.
  # This works because
  #   substitute(series, series_rk, t) = series[t] + lower order terms
  for o in 2:order
    for _t in RootedTreeIterator(o)
      t = copy(_t)
      series[t] += series_ex[t] - substitute(series, series_rk, t)
    end
  end

  return series
end


end # module
