module BSeriesSymPyExt

if isdefined(Base, :get_extension)
    using SymPy: SymPy
else
    using ..SymPy: SymPy
end

using BSeries: BSeries

function BSeries.compute_derivative(expression::SymPy.Sym{SymPy.PyCall.PyObject},
                                    variable::SymPy.Sym{SymPy.PyCall.PyObject})
    SymPy.diff(expression, variable)
end

function BSeries.latexify_default_dt(::Type{SymPy.Sym{SymPy.PyCall.PyObject}})
    SymPy.symbols("h", real = true)
end

end # module
