module BSeries

using Reexport: @reexport

@reexport using RootedTrees
using RootedTrees: RootedTree

@reexport using OrderedCollections: OrderedDict


export bseries, substitute

export modified_equation


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
function bseries(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)
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
    modified_equation(A, b, c, order)

Compute the B-series of the modified equation of the Runge-Kutta method with
Butcher coefficients `A, b, c` up to a prescribed `order`.
"""
function modified_equation(A, b, c, order)
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


end # module
