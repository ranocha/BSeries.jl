module BSeries

using Reexport: @reexport

@reexport using RootedTrees
using RootedTrees: RootedTree

@reexport using OrderedCollections: OrderedDict

@reexport using Symbolics: Symbolics, Num
using Symbolics: Differential, expand_derivatives


export bseries, substitute, compose

export modified_equation, modifying_integrator

export elementary_differentials


# This is a dirty workaround until the performance bugfix
# https://github.com/JuliaLang/julia/pull/42300
# is merged
@static if v"1.6" <= VERSION < v"1.7.0-rc2"
  Base.hastypemax(::Type{Bool}) = true
end


"""
    substitute(b, a, t::RootedTree)

Compute the coefficient corresponding to the tree `t` of the B-series that is
formed by substituting the B-series `b` into the B-series `a`.

# References

Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
substitute(b, a, t::RootedTree) = substitute!((nothing, nothing), b, a, t)

"""
    substitute!(cache, b, a, t::RootedTree)

Like [`substitute`](@ref) but with possible speed-up from memoization
(at the cost of more memory usage in the `cache`).
"""
function substitute!(cache, b, a, t::RootedTree)
  result = zero(first(values(a)) * first(values(b)))

  for (forest, skeleton) in PartitionIterator(t, cache)
    result += reduce(*, b[tree] for tree in forest) * a[skeleton]
  end

  return result
end


"""
    compose(b, a, t::RootedTree)

Compute the coefficient correspoding to the tree `t` of the B-series that is
formed by composing the B-series `a` with the B-series `b`.

# References

Section 3.1 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series
  Foundations of Computational Mathematics
  [DOI: 10.1007/s10208-010-9065-1](https://doi.org/10.1007/s10208-010-9065-1)
"""
function compose(b, a, t::RootedTree)
  result = zero(first(values(a)) * first(values(b)))

  for (forest, subtree) in SplittingIterator(t)
    if isempty(forest)
      result += a[subtree]
    else
      result += reduce(*, b[tree] for tree in forest) * a[subtree]
    end
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
  # Since the `keys` are ordered, we don't need to use nested loops of the form
  #   for o in 2:order
  #     for _t in RootedTreeIterator(o)
  #       t = copy(_t)
  # which are slightly less efficient due to additional computations and allocations.
  for t in keys(series)
    series[t] += series_rk[t] - substitute(series, series_ex, t)
  end

  return series
end


"""
    modified_equation(f, u, dt,
                      A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)

Compute the B-series of the modified equation of the Runge-Kutta method with
Butcher coefficients `A, b, c` up to the prescribed `order` with respect to
the ordinary differential equation ``u'(t) = f(u(t))`` with vector field `f`
and dependent variables `u` for a time step size `dt`.

Here, `u` is assumed to be a vector of variables from Symbolics.jl (`Symbolics.Num`)
and `f` is assumed to be a vector of expressions in these variables.

See Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series
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

  cache = create_cache(PartitionIterator, t)
  # cache = (first(cache), nothing) # TODO: Why is this faster?
  # cache = (nothing, nothing)

  # Recursively solve `substitute(series, series_rk, t) == series_ex[t]`.
  # This works because
  #   substitute(series, series_rk, t) = series[t] + lower order terms
  # Since the `keys` are ordered, we don't need to use nested loops of the form
  #   for o in 2:order
  #     for _t in RootedTreeIterator(o)
  #       t = copy(_t)
  # which are slightly less efficient due to additional computations and allocations.
  for t in keys(series)
    series[t] += series_ex[t] - substitute!(cache, series, series_rk, t)
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

Here, `u` is assumed to be a vector of variables from Symbolics.jl (`Symbolics.Num`)
and `f` is assumed to be a vector of expressions in these variables.

See Section 3.2 of
- Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)
  Algebraic Structures of B-series
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
    _compute_partial_derivatives!(d, f, u)
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

function _compute_partial_derivatives!(d, f, u)
  for idx in CartesianIndices(d)
    idx_tuple = Tuple(idx)
    f_idx = first(idx_tuple)
    u_idx = Base.tail(idx_tuple)

    # All this B-series analysis only really makes sense for smooth functions.
    # Hence, we can use the symmetry of the partial derivatives to speed-up
    # the computations - Hermann Amandus Schwarz helps us again!
    issorted(u_idx) || continue

    partial_derivative = f[f_idx]
    for i in u_idx
      partial_derivative = Differential(u[i])(partial_derivative)
    end
    d[idx] = expand_derivatives(partial_derivative)
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
