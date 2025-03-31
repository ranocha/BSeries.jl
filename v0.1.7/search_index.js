var documenterSearchIndex = {"docs":
[{"location":"benchmarks/#Benchmarks","page":"Benchmarks","title":"Benchmarks","text":"","category":"section"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"Here, we collect some simple benchmarks of BSeries.jl. Take them with a grain of salt since they run on virtual machines in the cloud to generate the documentation automatically. You can of course also copy the code and run the benchmarks locally yourself.","category":"page"},{"location":"benchmarks/#Comparing-different-symbolic-packages","page":"Benchmarks","title":"Comparing different symbolic packages","text":"","category":"section"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"Symbolic computations of modified_equations and modifying_integrators in BSeries.jl support","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"SymEngine.jl,\nSymPy.jl, and\nSymbolics.jl","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"as symbolic backends. Here, we compare them in the context of the explicit midpoint method and the nonlinear oscillator ODE","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"u(t) = frac1 u(t) ^2 beginpmatrix -u_2(t)  u_1(t) endpmatrix","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"This particular combination of explicit Runge-Kutta method and ODE is special since the explicit midpoint method is unconditionally energy-conserving for this problem[RanochaKetcheson2020].","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"First, we set up some code to perform the benchmarks. Here, we use a very naive approach, run the code twice (to see the effect of compilation) and use @time to print the runtime. More sophisticated approaches should of course use something like @benchmark from BenchmarkTools.jl. However, this simple and cheap version suffices to compare the orders of magnitude.","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"using BSeries, StaticArrays\n\nfunction benchmark(u, dt, subs, order)\n  # explicit midpoint method\n  A = @SArray [0 0; 1//2 0]\n  b = @SArray [0, 1//1]\n  c = @SArray [0, 1//2]\n\n  # nonlinear oscillator\n  f = [-u[2], u[1]] / (u[1]^2 + u[2]^2)\n\n  println(\"\\n Computing the series coefficients:\")\n  @time coefficients = modifying_integrator(A, b, c, order)\n  @time coefficients = modifying_integrator(A, b, c, order)\n\n  println(\"\\n Computing the series including elementary differentials:\")\n  @time series = modifying_integrator(f, u, dt, A, b, c, order)\n  @time series = modifying_integrator(f, u, dt, A, b, c, order)\n\n  substitution_variables = Dict(u[1] => 1//1, u[2] => 0//1)\n\n  println(\"\\n Substituting the initial condition:\")\n  @time subs.(series, (substitution_variables, ))\n  @time subs.(series, (substitution_variables, ))\n\n  println(\"\\n\")\nend","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"Next, we load the symbolic packages and run the benchmarks.","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"using SymPy # generates annoying output online when conda installs sympy","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"using SymEngine: SymEngine\nusing SymPy: SymPy\nusing Symbolics: Symbolics\n\nprintln(\"SymEngine\")\ndt   = SymEngine.symbols(\"dt\")\nu    = SymEngine.symbols(\"u1, u2\")\nsubs = SymEngine.subs\nbenchmark(u, dt, subs, 8)\n\nprintln(\"SymPy\")\ndt   = SymPy.symbols(\"dt\")\nu    = SymPy.symbols(\"u1, u2\")\nsubs = SymPy.subs\nbenchmark(u, dt, subs, 8)\n\nprintln(\"Symbolics\")\nSymbolics.@variables dt\nu = Symbolics.@variables u1 u2\nsubs = Symbolics.substitute\nbenchmark(u, dt, subs, 8)","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"These results were obtained using the following versions.","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"using InteractiveUtils\nversioninfo()\n\nusing Pkg\nPkg.status([\"SymEngine\", \"SymPy\", \"Symbolics\"])\nnothing # hide","category":"page"},{"location":"benchmarks/#References","page":"Benchmarks","title":"References","text":"","category":"section"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"[RanochaKetcheson2020]: ","category":"page"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"Hendrik Ranocha and David Ketcheson (2020)   Energy Stability of Explicit Runge-Kutta Methods for Nonautonomous or   Nonlinear Problems.   SIAM Journal on Numerical Analysis   DOI: 10.1137/19M1290346","category":"page"},{"location":"api_reference/#BSeries.jl-API","page":"API reference","title":"BSeries.jl API","text":"","category":"section"},{"location":"api_reference/","page":"API reference","title":"API reference","text":"CurrentModule = BSeries","category":"page"},{"location":"api_reference/","page":"API reference","title":"API reference","text":"Modules = [BSeries]","category":"page"},{"location":"api_reference/#BSeries.ExactSolution","page":"API reference","title":"BSeries.ExactSolution","text":"ExactSolution{T}()\n\nLazy representation of the B-series of the exact solution of an ordinary differential equation using coefficients of type at least as representative as T.\n\n\n\n\n\n","category":"type"},{"location":"api_reference/#BSeries.bseries-Tuple{AbstractMatrix{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, Any}","page":"API reference","title":"BSeries.bseries","text":"bseries(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)\n\nCompute the B-series of the  Runge-Kutta method with Butcher coefficients A, b, c up to a prescribed order\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.compose-Tuple{Any, Any, RootedTrees.RootedTree}","page":"API reference","title":"BSeries.compose","text":"compose(b, a, t::RootedTree)\n\nCompute the coefficient correspoding to the tree t of the B-series that is formed by composing the B-series a with the B-series b.\n\nReferences\n\nSection 3.1 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.compute_derivative","page":"API reference","title":"BSeries.compute_derivative","text":"compute_derivative(expression, variable)\n\nInternal function specialized on symbolic variables and expressions from\n\nSymEngine.jl,\nSymPy.jl, and\nSymbolics.jl\n\nif these packages are loaded (via Requires.jl).\n\n\n\n\n\n","category":"function"},{"location":"api_reference/#BSeries.elementary_differentials-Tuple{Any, Any, Any}","page":"API reference","title":"BSeries.elementary_differentials","text":"elementary_differentials(f, u, order)\n\nCompute all elementary differentials of the vector field f with independent variables u up to the given order. The return value can be indexed by rooted trees to obtain the corresponding elementary differential.\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.modified_equation-Tuple{AbstractMatrix{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, Any}","page":"API reference","title":"BSeries.modified_equation","text":"modified_equation(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)\n\nCompute the B-series of the modified equation of the Runge-Kutta method with Butcher coefficients A, b, c up to the prescribed order.\n\nGiven an ordinary differential equation (ODE) u(t) = f(u(t)) and a Runge-Kutta method, the idea is to interpret the numerical solution with given time step size as exact solution of a modified ODE u(t) = fₕ(u(t)).\n\nnote: Normalization by elementary differentials\nThe coefficients of the B-series returned by this method need to be multiplied by the corresponding elementary differential of the input vector field f.\n\nReferences\n\nSection 3.2 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.modified_equation-Tuple{Any, Any, Any, AbstractMatrix{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, Any}","page":"API reference","title":"BSeries.modified_equation","text":"modified_equation(f, u, dt,\n                  A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)\n\nCompute the B-series of the modified equation of the Runge-Kutta method with Butcher coefficients A, b, c up to the prescribed order with respect to the ordinary differential equation u(t) = f(u(t)) with vector field f and dependent variables u for a time step size dt.\n\nHere, u is assumed to be a vector of symbolic variables and f is assumed to be a vector of expressions in these variables. Currently, symbolic variables from\n\nSymEngine.jl,\nSymPy.jl, and\nSymbolics.jl\n\nare supported.\n\nReferences\n\nSection 3.2 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.modifying_integrator-Tuple{AbstractMatrix{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, Any}","page":"API reference","title":"BSeries.modifying_integrator","text":"modifying_integrator(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)\n\nCompute the B-series of a \"modifying integrator\" equation of the Runge-Kutta method with Butcher coefficients A, b, c up to the prescribed order.\n\nGiven an ordinary differential equation (ODE) u(t) = f(u(t)) and a Runge-Kutta method, the idea is to find a modified ODE u(t) = fₕ(u(t)) such that the numerical solution with given time step size is the exact solution of the original ODE.\n\nnote: Normalization by elementary differentials\nThe coefficients of the B-series returned by this method need to be multiplied by the corresponding elementary differential of the input vector field f.\n\nReferences\n\nSection 3.2 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.modifying_integrator-Tuple{Any, Any, Any, AbstractMatrix{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, Any}","page":"API reference","title":"BSeries.modifying_integrator","text":"modifying_integrator(f, u, dt,\n                     A::AbstractMatrix, b::AbstractVector, c::AbstractVector, order)\n\nCompute the B-series of a \"modifying integrator\" equation of the Runge-Kutta method with Butcher coefficients A, b, c up to the prescribed order with respect to the ordinary differential equation u(t) = f(u(t)) with vector field f and dependent variables u for a time step size dt.\n\nHere, u is assumed to be a vector of symbolic variables and f is assumed to be a vector of expressions in these variables. Currently, symbolic variables from\n\nSymEngine.jl,\nSymPy.jl, and\nSymbolics.jl\n\nare supported.\n\nReferences\n\nSection 3.2 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"api_reference/#BSeries.substitute-Tuple{Any, Any, RootedTrees.RootedTree}","page":"API reference","title":"BSeries.substitute","text":"substitute(b, a, t::RootedTree)\n\nCompute the coefficient correspoding to the tree t of the B-series that is formed by substituting the B-series b into the B-series a.\n\nReferences\n\nSection 3.2 of\n\nPhilippe Chartier, Ernst Hairer, Gilles Vilmart (2010) Algebraic Structures of B-series. Foundations of Computational Mathematics DOI: 10.1007/s10208-010-9065-1\n\n\n\n\n\n","category":"method"},{"location":"license/","page":"License","title":"License","text":"EditURL = \"https://github.com/ranocha/BSeries.jl/blob/main/LICENSE.md\"","category":"page"},{"location":"license/#License","page":"License","title":"License","text":"","category":"section"},{"location":"license/","page":"License","title":"License","text":"MIT LicenseCopyright (c) 2021-present Hendrik RanochaPermission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.","category":"page"},{"location":"tutorials/modifying_integrators/#Modifying-integrators","page":"Modifying integrators","title":"Modifying integrators","text":"","category":"section"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"This tutorial describes the API of BSeries.jl related to the notion of modifying integrators. The main API entry point is modifying_integrator.","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"Given a first-order autonomous ordinary differential equation (ODE)","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"u(t) = f(u(t))","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"and a B-series time integration method, the idea is to find a modified ODE","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"u(t) = f_h(u(t))","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"such that the numerical solution with given time step size h of the original ODE is the exact solution of the modified ODE, see [ChartierHairerVilmart2007] and [ChartierHairerVilmart2010].","category":"page"},{"location":"tutorials/modifying_integrators/#Lotka-Volterra-model","page":"Modifying integrators","title":"Lotka-Volterra model","text":"","category":"section"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"Here, we consider the explicit Euler method to solve the classical Lotka-Volterra model","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"p(t) = (2 - q) p\nquad\nq(t) = (p - 1) q","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"First, we set up the ODE and compute some numerical solutions using OrdinaryDiffEq.jl.","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"using OrdinaryDiffEq, StaticArrays\n\nfunction f(u, params, t)\n  p, q = u\n  dp = (2 - q) * p\n  dq = (p - 1) * q\n  return SVector(dp, dq)\nend\n\nu0 = SVector(1.5, 2.25)\ntspan = (0.0, 15.0)\node = ODEProblem(f, u0, tspan)\n\ndt = 0.35\nsol_euler = solve(ode, Euler(), dt=dt)\nsol_ref = solve(ode, Tsit5())\nnothing # hide","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"Next, we look at some phase space plots of the numerical solution.","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"using LaTeXStrings, Plots\n\nfig = plot(xguide=L\"$q$\", yguide=L\"$p$\")\ndefault(linewidth=2)\nplot!(fig, sol_ref, vars=(2, 1), label=\"Reference solution\")\nscatter!(fig, last.(sol_euler.u), first.(sol_euler.u),\n         label=\"Explicit Euler, dt = $dt\")\nplot!(fig, xlims=(0.0, 4.0), ylims=(0.0, 2.5))\n\nsavefig(fig, \"lotka_volterra_original.svg\"); nothing # hide","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"(Image: )","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"The exact solution of this problem is periodic, but the explicit Euler method produces an unstable trajectory. Here, we used an especially large time step to more clearly illustrate what will follow, but the qualitative behavior is the same for any time step size.","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"Next, we will derive a \"modifying integrator\". What this means is that we will determine a perturbed ODE right-hand side (RHS) such that when Euler's method is applied to the perturbed RHS, the result is the exact solution to the original Lotka-Volterra system. The perturbed system takes the form of a power series in the time step size dt, and in order to compute with it we will truncate it at a certain order. We can compare the accuracy (and qualitative behavior) obtained by truncating at different orders.","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"Here, we use Symbolics.jl for the symbolic computations.","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"using BSeries, Symbolics\n\n# Explicit Euler method\nA = @SArray [0//1;]\nb = @SArray [1//1]\nc = @SArray [0//1]\n\n# Setup of symbolic variables\n@variables dt_sym\nu_sym = @variables p q\nf_sym = [f(u_sym, nothing, nothing)...]\n\nfor truncation_order in 2:4\n  series = modifying_integrator(f_sym, u_sym, dt_sym, A, b, c, truncation_order)\n  series = Symbolics.substitute.(series, dt_sym => dt)\n  modified_f, _ = build_function(series, u_sym, expression=Val(false))\n  modified_ode = ODEProblem((u, params, t) -> modified_f(u), ode.u0, tspan)\n  modified_sol_euler = solve(modified_ode, Euler(), dt=dt)\n  plot!(fig, modified_sol_euler, vars=(2, 1),\n        label=\"Euler, modified ODE order $(truncation_order-1)\")\nend\nplot!(fig, xlims=(0.0, 4.0), ylims=(0.0, 2.5))\n\nsavefig(fig, \"lotka_volterra_modified.svg\"); nothing # hide","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"(Image: )","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"We see that if we include one additional term, the resulting trajectory still grows, while with two additional terms the solution appears to be dissipative. With each additional term, the solution gets closer to the exact solution of the original problem, and with three added terms it is hard to see the difference between them at this scale.","category":"page"},{"location":"tutorials/modifying_integrators/#References","page":"Modifying integrators","title":"References","text":"","category":"section"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"[ChartierHairerVilmart2007]: ","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"Philippe Chartier, Ernst Hairer and Gilles Vilmart (2007)   Numerical integrators based on modified differential equations   DOI: 10.1090/S0025-5718-07-01967-9","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"[ChartierHairerVilmart2010]: ","category":"page"},{"location":"tutorials/modifying_integrators/","page":"Modifying integrators","title":"Modifying integrators","text":"Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)   Algebraic Structures of B-series.   Foundations of Computational Mathematics   DOI: 10.1007/s10208-010-9065-1","category":"page"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"EditURL = \"https://github.com/ranocha/BSeries.jl/blob/main/CONTRIBUTING.md\"","category":"page"},{"location":"contributing/#Contributing","page":"Contributing","title":"Contributing","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"ContributingBSeries.jl is an open-source project and we are very happy to accept contributions from the community. Please feel free to open issues or submit patches (preferably as pull requests) any time. For planned larger contributions, it is often beneficial to get in contact first, for example via issues.BSeries.jl and its contributions are licensed under the MIT license (see License). As a contributor, you certify that all your contributions are in conformance with the Developer Certificate of Origin (Version 1.1), which is reproduced below.Developer Certificate of Origin (Version 1.1)The following text was taken from https://developercertificate.org:Developer Certificate of Origin\nVersion 1.1\n\nCopyright (C) 2004, 2006 The Linux Foundation and its contributors.\n1 Letterman Drive\nSuite D4700\nSan Francisco, CA, 94129\n\nEveryone is permitted to copy and distribute verbatim copies of this\nlicense document, but changing it is not allowed.\n\n\nDeveloper's Certificate of Origin 1.1\n\nBy making a contribution to this project, I certify that:\n\n(a) The contribution was created in whole or in part by me and I\n    have the right to submit it under the open source license\n    indicated in the file; or\n\n(b) The contribution is based upon previous work that, to the best\n    of my knowledge, is covered under an appropriate open source\n    license and I have the right under that license to submit that\n    work with modifications, whether created in whole or in part\n    by me, under the same open source license (unless I am\n    permitted to submit under a different license), as indicated\n    in the file; or\n\n(c) The contribution was provided directly to me by some other\n    person who certified (a), (b) or (c) and I have not modified\n    it.\n\n(d) I understand and agree that this project and the contribution\n    are public and that a record of the contribution (including all\n    personal information I submit with it, including my sign-off) is\n    maintained indefinitely and may be redistributed consistent with\n    this project or the open source license(s) involved.","category":"page"},{"location":"#BSeries.jl","page":"Home","title":"BSeries.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The Julia library BSeries.jl is work in progress. Nevertheless, we follow semantic versioning. Thus, you can safely use the package right now. Extended documentation will be provided in the future.","category":"page"},{"location":"","page":"Home","title":"Home","text":"BSeries.jl re-exports everything from RootedTrees.jl. However, if you rely on functionality from that package, you should also include it explicitly in your project dependencies to track breaking changes, since the version numbers of RootedTrees.jl and BSeries.jl are not necessarily synchronized.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"BSeries.jl is a registered Julia package. Thus, you can install it from the Julia REPL via","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.add(\"BSeries\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you want to update BSeries.jl, you can use","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.update(\"BSeries\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"As usual, if you want to update BSeries.jl and all other packages in your current project, you can execute","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.update()","category":"page"},{"location":"#Referencing","page":"Home","title":"Referencing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you use BSeries.jl for your research, please cite it using the bibtex entry","category":"page"},{"location":"","page":"Home","title":"Home","text":"@misc{ranocha2021bseries,\n  title={{BSeries.jl}: {C}omputing with {B}-series in {J}ulia},\n  author={Ranocha, Hendrik and Ketcheson, David I},\n  year={2021},\n  month={09},\n  howpublished={\\url{https://github.com/ranocha/BSeries.jl}},\n  doi={10.5281/zenodo.5534602}\n}","category":"page"},{"location":"","page":"Home","title":"Home","text":"Please also cite the appropriate references for specific functions you use, which can be obtained from their docstrings.","category":"page"},{"location":"#License-and-contributing","page":"Home","title":"License and contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This project is licensed under the MIT license (see License). Since it is an open-source project, we are very happy to accept contributions from the community. Please refer to the section Contributing for more details.","category":"page"},{"location":"tutorials/modified_equations/#Modified-equations","page":"Modified equations","title":"Modified equations","text":"","category":"section"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"This tutorial describes the API of BSeries.jl related to the notion of modified equations. The main API entry point is modified_equation.","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"Given a first-order autonomous ordinary differential equation (ODE)","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"u(t) = f(u(t))","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"and a B-series time integration method, the idea is to interpret the numerical solution with given time step size h of the original ODE as the exact solution of the modified ODE","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"u(t) = f_h(u(t))","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"see for example [ChartierHairerVilmart2010].","category":"page"},{"location":"tutorials/modified_equations/#Lotka-Volterra-model","page":"Modified equations","title":"Lotka-Volterra model","text":"","category":"section"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"Here, we reproduce the example on p. 340 of [HairerLubichWanner2006]. Thus, we consider the explicit Euler method to solve the classical Lotka-Volterra model","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"p(t) = (2 - q) p\nquad\nq(t) = (p - 1) q","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"First, we set up the ODE and compute some numerical solutions using OrdinaryDiffEq.jl.","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"using OrdinaryDiffEq, StaticArrays\n\nfunction f(u, params, t)\n  p, q = u\n  dp = (2 - q) * p\n  dq = (p - 1) * q\n  return SVector(dp, dq)\nend\n\nu0 = SVector(1.5, 2.25)\ntspan = (0.0, 15.0)\node = ODEProblem(f, u0, tspan)\n\ndt = 0.1\nsol_euler = solve(ode, Euler(), dt=dt)\nsol_ref = solve(ode, Tsit5())\nnothing # hide","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"Next, we look at some phase space plots of the numerical solution.","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"using LaTeXStrings, Plots\n\nfig = plot(xguide=L\"$q$\", yguide=L\"$p$\")\ndefault(linewidth=2)\nplot!(fig, sol_ref, vars=(2, 1), label=\"Reference solution\")\nscatter!(fig, last.(sol_euler.u), first.(sol_euler.u),\n         label=\"Explicit Euler, dt = $dt\")\nplot!(fig, xlims=(0.0, 9.0), ylims=(0.0, 5.0))\n\nsavefig(fig, \"lotka_volterra_original.svg\"); nothing # hide","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"(Image: )","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"The exact solution of this problem is periodic, but the explicit Euler method produces an unstable trajectory. Here, we used an especially large time step to more clearly illustrate what will follow, but the qualitative behavior is the same for any time step size.","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"Next, we will derive the \"modified equation\" of the explicit Euler method and solve this new ODE to high accuracy. The perturbed system takes the form of a power series in the time step size dt, and in order to compute with it we will truncate it at a certain order.","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"Here, we use Symbolics.jl for the symbolic computations.","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"using BSeries, Symbolics\n\nfunction solve_modified_equation(ode, truncation_orders, dt)\n  # Explicit Euler method\n  A = @SArray [0//1;]\n  b = @SArray [1//1]\n  c = @SArray [0//1]\n\n  # Setup of symbolic variables\n  @variables dt_sym\n  u_sym = @variables p q\n  f_sym = [f(u_sym, nothing, nothing)...]\n\n  sol_euler = solve(ode, Euler(), dt=dt)\n\n  fig = plot(xguide=L\"$q$\", yguide=L\"$p$\")\n  default(linewidth=2)\n  scatter!(fig, last.(sol_euler.u), first.(sol_euler.u),\n          label=\"Explicit Euler, dt = $dt\")\n\n  for truncation_order in truncation_orders\n    series = modified_equation(f_sym, u_sym, dt_sym, A, b, c, truncation_order)\n    series = Symbolics.substitute.(series, dt_sym => dt)\n    modified_f, _ = build_function(series, u_sym, expression=Val(false))\n    modified_ode = ODEProblem((u, params, t) -> modified_f(u), ode.u0, tspan)\n    modified_sol = solve(modified_ode, Tsit5())\n    plot!(fig, modified_sol, vars=(2, 1),\n          label=\"Modified ODE, order $(truncation_order-1)\")\n  end\n  fig\nend\n\nfig = solve_modified_equation(ode, 2, dt)\n\nsavefig(fig, \"lotka_volterra_modified1.svg\"); nothing # hide","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"(Image: )","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"The exact solution of the Lotka-Volterra model is periodic, but Euler's method generates a solution with growing amplitude. The modified equations accurately predict this.","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"Now we go to the next order and increase the time step size dt slightly.","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"fig = solve_modified_equation(ode, 2:3, 0.11)\n\nsavefig(fig, \"lotka_volterra_modified2.svg\"); nothing # hide","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"(Image: )","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"Using a larger step size, we see that the first-order modified equations are not fully accurate, but by including the O(h^2) terms we get much better accuracy at late times. Let's keep going.","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"fig = solve_modified_equation(ode, 2:4, 0.12)\n\nsavefig(fig, \"lotka_volterra_modified3.svg\"); nothing # hide","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"(Image: )","category":"page"},{"location":"tutorials/modified_equations/#References","page":"Modified equations","title":"References","text":"","category":"section"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"[HairerLubichWanner2006]: ","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"Ernst Hairer, Christian Lubich, Gerhard Wanner (2006)   Geometric Numerical Integration.   DOI: 10.1007/3-540-30666-8","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"[ChartierHairerVilmart2010]: ","category":"page"},{"location":"tutorials/modified_equations/","page":"Modified equations","title":"Modified equations","text":"Philippe Chartier, Ernst Hairer, Gilles Vilmart (2010)   Algebraic Structures of B-series.   Foundations of Computational Mathematics   DOI: 10.1007/s10208-010-9065-1","category":"page"}]
}
