module SymbolicsExt

if isdefined(Base, :get_extension)
    using Symbolics: Symbolics
else
    using ..Symbolics: Symbolics
end

using BSeries: BSeries

function BSeries.compute_derivative(expression::Symbolics.Num, variable::Symbolics.Num)
    Symbolics.expand_derivatives(Symbolics.Differential(variable)(expression))
end

BSeries.latexify_default_dt(::Type{Symbolics.Num}) = (Symbolics.@variables(h); h)

end # module
