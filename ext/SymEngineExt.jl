module SymEngineExt

if isdefined(Base, :get_extension)
    using SymEngine: SymEngine
else
    using ..SymEngine: SymEngine
end

using BSeries: BSeries

function BSeries.compute_derivative(expression::SymEngine.Basic, variable::SymEngine.Basic)
    SymEngine.diff(expression, variable)
end

BSeries.latexify_default_dt(::Type{SymEngine.Basic}) = SymEngine.symbols("h")

end # module
