module BSeries

using Reexport: @reexport

@reexport using RootedTrees
using RootedTrees: RootedTree

@reexport using OrderedCollections: OrderedDict


export substitute


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


end # module
