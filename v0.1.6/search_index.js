var documenterSearchIndex = {"docs":
[{"location":"api_reference/#BSeries.jl-API","page":"API reference","title":"BSeries.jl API","text":"","category":"section"},{"location":"api_reference/","page":"API reference","title":"API reference","text":"CurrentModule = BSeries","category":"page"},{"location":"api_reference/","page":"API reference","title":"API reference","text":"Modules = [BSeries]","category":"page"},{"location":"api_reference/#BSeries.ExactSolution","page":"API reference","title":"BSeries.ExactSolution","text":"ExactSolution{T}()\n\nLazy representation of the B-series of the exact solution of an ordinary differential equation using coefficients of type at least as representative as T.\n\n\n\n\n\n","category":"type"},{"location":"api_reference/#BSeries.bseries-Tuple{AbstractMatrix{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, Any}","page":"API reference","title":"BSeries.bseries","text":"bseries(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)\n\nCompute the B-series of the  Runge-Kutta method with Butcher coefficients A, b, c up to a prescribed order\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.compose-Tuple{Any, Any, RootedTrees.RootedTree}","page":"API reference","title":"BSeries.compose","text":"compose(b, a, t::RootedTree)\n\nCompute the coefficient correspoding to the tree t of the B-series that is formed by composing the B-series a with the B-series b.\n\nReferences\n\nSection 3.1 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.compute_derivative","page":"API reference","title":"BSeries.compute_derivative","text":"compute_derivative(expression, variable)\n\nInternal function specialized on symbolic variables and expressions from\n\nSymEngine.jl,\nSymPy.jl, and\nSymbolics.jl\n\nif these packages are loaded (via Requires.jl).\n\n\n\n\n\n","category":"function"},{"location":"api_reference/#BSeries.elementary_differentials-Tuple{Any, Any, Any}","page":"API reference","title":"BSeries.elementary_differentials","text":"elementary_differentials(f, u, order)\n\nCompute all elementary differentials of the vector field f with independent variables u up to the given order. The return value can be indexed by rooted trees to obtain the corresponding elementary differential.\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.modified_equation-Tuple{AbstractMatrix{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, Any}","page":"API reference","title":"BSeries.modified_equation","text":"modified_equation(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)\n\nCompute the B-series of the modified equation of the Runge-Kutta method with Butcher coefficients A, b, c up to the prescribed order.\n\nGiven an ordinary differential equation (ODE) u(t) = f(u(t)) and a Runge-Kutta method, the idea is to interpret the numerical solution with given time step size as exact solution of a modified ODE u(t) = fₕ(u(t)).\n\nnote: Normalization by elementary differentials\nThe coefficients of the B-series returned by this method need to be multiplied by the corresponding elementary differential of the input vector field f.\n\nReferences\n\nSection 3.2 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.modified_equation-Tuple{Any, Any, Any, AbstractMatrix{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, Any}","page":"API reference","title":"BSeries.modified_equation","text":"modified_equation(f, u, dt,\n                  A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)\n\nCompute the B-series of the modified equation of the Runge-Kutta method with Butcher coefficients A, b, c up to the prescribed order with respect to the ordinary differential equation u(t) = f(u(t)) with vector field f and dependent variables u for a time step size dt.\n\nHere, u is assumed to be a vector of symbolic variables and f is assumed to be a vector of expressions in these variables. Currently, symbolic variables from\n\nSymEngine.jl,\nSymPy.jl, and\nSymbolics.jl\n\nare supported.\n\nReferences\n\nSection 3.2 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.modifying_integrator-Tuple{AbstractMatrix{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, Any}","page":"API reference","title":"BSeries.modifying_integrator","text":"modifying_integrator(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)\n\nCompute the B-series of a \"modifying integrator\" equation of the Runge-Kutta method with Butcher coefficients A, b, c up to the prescribed order.\n\nGiven an ordinary differential equation (ODE) u(t) = f(u(t)) and a Runge-Kutta method, the idea is to find a modified ODE u(t) = fₕ(u(t)) such that the numerical solution with given time step size is the exact solution of the original ODE.\n\nnote: Normalization by elementary differentials\nThe coefficients of the B-series returned by this method need to be multiplied by the corresponding elementary differential of the input vector field f.\n\nReferences\n\nSection 3.2 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.modifying_integrator-Tuple{Any, Any, Any, AbstractMatrix{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, Any}","page":"API reference","title":"BSeries.modifying_integrator","text":"modifying_integrator(f, u, dt,\n                     A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)\n\nCompute the B-series of a \"modifying integrator\" equation of the Runge-Kutta method with Butcher coefficients A, b, c up to the prescribed order with respect to the ordinary differential equation u(t) = f(u(t)) with vector field f and dependent variables u for a time step size dt.\n\nHere, u is assumed to be a vector of symbolic variables and f is assumed to be a vector of expressions in these variables. Currently, symbolic variables from\n\nSymEngine.jl,\nSymPy.jl, and\nSymbolics.jl\n\nare supported.\n\nReferences\n\nSection 3.2 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.substitute-Tuple{Any, Any, RootedTrees.RootedTree}","page":"API reference","title":"BSeries.substitute","text":"substitute(b, a, t::RootedTree)\n\nCompute the coefficient correspoding to the tree t of the B-series that is formed by substituting the B-series b into the B-series a.\n\nReferences\n\nSection 3.2 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"license/","page":"License","title":"License","text":"EditURL = \"https://github.com/ranocha/BSeries.jl/blob/main/LICENSE.md\"","category":"page"},{"location":"license/#License","page":"License","title":"License","text":"","category":"section"},{"location":"license/","page":"License","title":"License","text":"MIT LicenseCopyright (c) 2021-present Hendrik RanochaPermission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.","category":"page"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"EditURL = \"https://github.com/ranocha/BSeries.jl/blob/main/CONTRIBUTING.md\"","category":"page"},{"location":"contributing/#Contributing","page":"Contributing","title":"Contributing","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"ContributingBSeries.jl is an open-source project and we are very happy to accept contributions from the community. Please feel free to open issues or submit patches (preferably as pull requests) any time. For planned larger contributions, it is often beneficial to get in contact first, for example via issues.BSeries.jl and its contributions are licensed under the MIT license (see License). As a contributor, you certify that all your contributions are in conformance with the Developer Certificate of Origin (Version 1.1), which is reproduced below.Developer Certificate of Origin (Version 1.1)The following text was taken from https://developercertificate.org:Developer Certificate of Origin\nVersion 1.1\n\nCopyright (C) 2004, 2006 The Linux Foundation and its contributors.\n1 Letterman Drive\nSuite D4700\nSan Francisco, CA, 94129\n\nEveryone is permitted to copy and distribute verbatim copies of this\nlicense document, but changing it is not allowed.\n\n\nDeveloper's Certificate of Origin 1.1\n\nBy making a contribution to this project, I certify that:\n\n(a) The contribution was created in whole or in part by me and I\n    have the right to submit it under the open source license\n    indicated in the file; or\n\n(b) The contribution is based upon previous work that, to the best\n    of my knowledge, is covered under an appropriate open source\n    license and I have the right under that license to submit that\n    work with modifications, whether created in whole or in part\n    by me, under the same open source license (unless I am\n    permitted to submit under a different license), as indicated\n    in the file; or\n\n(c) The contribution was provided directly to me by some other\n    person who certified (a), (b) or (c) and I have not modified\n    it.\n\n(d) I understand and agree that this project and the contribution\n    are public and that a record of the contribution (including all\n    personal information I submit with it, including my sign-off) is\n    maintained indefinitely and may be redistributed consistent with\n    this project or the open source license(s) involved.","category":"page"},{"location":"#BSeries.jl","page":"Home","title":"BSeries.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The Julia library BSeries.jl is work in progress. Nevertheless, we follow semantic versioning. Thus, you can safely use the package right now. Extended documentation will be provided in the future.","category":"page"},{"location":"","page":"Home","title":"Home","text":"BSeries.jl re-exports everything from RootedTrees.jl. However, if you rely on functionality from that package, you should also include it explicitly in your project dependencies to track breaking changes, since the version numbers of RootedTrees.jl and BSeries.jl are not necessarily synchronized.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"BSeries.jl is a registered Julia package. Thus, you can install it from the Julia REPL via","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.add(\"BSeries\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you want to update BSeries.jl, you can use","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.update(\"BSeries\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"As usual, if you want to update BSeries.jl and all other packages in your current project, you can execute","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.update()","category":"page"},{"location":"#Referencing","page":"Home","title":"Referencing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you use BSeries.jl for your research, please cite it using the bibtex entry","category":"page"},{"location":"","page":"Home","title":"Home","text":"@misc{ranocha2021bseries,\n  title={{BSeries.jl}: {C}omputing with {B}-series in {J}ulia},\n  author={Ranocha, Hendrik and Ketcheson, David I},\n  year={2021},\n  month={09},\n  howpublished={\\url{https://github.com/ranocha/BSeries.jl}},\n  doi={10.5281/zenodo.5534602}\n}","category":"page"},{"location":"","page":"Home","title":"Home","text":"Please also cite the appropriate references for specific functions you use, which can be obtained from their docstrings.","category":"page"},{"location":"#License-and-contributing","page":"Home","title":"License and contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This project is licensed under the MIT license (see License). Since it is an open-source project, we are very happy to accept contributions from the community. Please refer to the section Contributing for more details.","category":"page"}]
}
