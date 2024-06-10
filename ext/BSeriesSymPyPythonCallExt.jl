module BSeriesSymPyPythonCallExt

if isdefined(Base, :get_extension)
    using SymPyPythonCall: SymPyPythonCall
else
    using ..SymPyPythonCall: SymPyPythonCall
end

using BSeries: BSeries

function BSeries.compute_derivative(expression::SymPyPythonCall.Sym,
                                    variable::SymPyPythonCall.Sym)
    SymPyPythonCall.diff(expression, variable)
end

function BSeries.latexify_default_dt(::Type{SymPyPythonCall.Sym})
    SymPyPythonCall.symbols("h", real = true)
end

end # module
