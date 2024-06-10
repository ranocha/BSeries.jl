module BSeriesSymPyExt

if isdefined(Base, :get_extension)
    using SymPy: SymPy
else
    using ..SymPy: SymPy
end

using BSeries: BSeries

function BSeries.compute_derivative(expression::SymPy.Sym, variable::SymPy.Sym)
    SymPy.diff(expression, variable)
end

BSeries.latexify_default_dt(::Type{SymPy.Sym}) = SymPy.symbols("h", real = true)

end # module
