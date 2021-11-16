
struct LatexifyElementaryDifferential{T<:RootedTree}
  t::T
  f::String
end

Latexify.@latexrecipe function _(led::LatexifyElementaryDifferential)
  return LaTeXString("F_{" * led.f * "}\\mathopen{}\\left( " *
                     Latexify.latexify(led.t) * " \\right)\\mathclose{}")
end

latexify_default_dt(type) = "h"

Latexify.@latexrecipe function _(series::TruncatedBSeries;
                                 f="f", dt="no sensible default value",
                                 reduce_order_by=0)

  # `@latexrecipe` turns many assignments to keyword arguments into expressions,
  # which doesn't work at all. Hence, we need to use a recognized type and
  # set sensible default values here.
  if dt == "no sensible default value"
    dt = latexify_default_dt(valtype(series))
  end

  expressions = []
  for (t, val) in series
    iszero(val) && continue
    elementary_differential = LatexifyElementaryDifferential(t, f)
    if dt isa Symbol || dt isa AbstractString
      # insert the symbol of dt
      val_symmetry = val / symmetry(t)
      reduced_order = order(t) - reduce_order_by
      if isone(val_symmetry)
        if iszero(reduced_order)
          push!(expressions,
                :($(elementary_differential)))
        elseif isone(reduced_order)
          push!(expressions,
                :($dt *
                  $(elementary_differential)))
        else
          push!(expressions,
                :($dt^$(reduced_order) *
                  $(elementary_differential)))
        end
      else
        push!(expressions,
              :($(val / symmetry(t)) *
                $dt^$(order(t) - reduce_order_by) *
                $(elementary_differential)))
      end
    else
      # assume that we can do arithmetic with dt
      factor = val / symmetry(t) * dt^(order(t) - reduce_order_by)
      iszero(factor) && continue
      if isone(factor)
        push!(expressions,
              :($(elementary_differential)))
      else
        push!(expressions,
              :($(factor) *
                $(elementary_differential)))
      end
    end
  end
  result = expressions[1]
  for i in 2:length(expressions)
    result = :($result + $(expressions[i]))
  end

  return result
end
