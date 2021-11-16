
latexify_default_dt(::Type{Symbolics.Num}) = (Symbolics.@variables(h); h)
