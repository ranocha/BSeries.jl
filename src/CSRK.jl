#Packages which would be useful eventually

using BSeries
using SymPy; 

using LinearAlgebra
using PyCall



M = [1 2 0 1 1;
     4 5 6 4 5;
     7 8 9 0 1;
     1 3 4 3 1;
     1 0 2 1 0]

M = [ 1 2 3;
     2 3 4;
     1 0 2]

"""
    elementary_differentials_csrk(M,tree)

This function calculates the elementary differential for a given 
square matrix 'M' and a given RootedTree.

"""
function elementary_differentials_csrk(M,rootedtree)
    tree = rootedtree.level_sequence
    m = maximum(tree)
    l = length(tree)
    #create the variables
    variables = []
    for i in 1:m
        var_name = "x$i"
        var = Sym(var_name)
        push!(variables, var)
    end
    inverse_counter = l-1
    #stablish initial integrand, which is the rightmost leaf (last node)
    if l > 1
        integrand = integrate(PolinomialA(M,variables[tree[end]-1],variables[tree[end]]),(variables[tree[end]],0,1))
    else
        return 1
    end
    while inverse_counter > 1 
        #println("A_",variables[tree[inverse_counter]-1],variables[tree[inverse_counter]])
        pseudo_integrand = PolinomialA(M,variables[tree[inverse_counter]-1],variables[tree[inverse_counter]])*integrand
        integrand = integrate(pseudo_integrand,(variables[tree[inverse_counter]],0,1))
        #println(integrand)
        inverse_counter -= 1
    end
    #println("B_x1")
    #multiply for the Basis_Polonimial, i.e. the Polinomial B
    #return the integral with respect to x1.
    return integrate(PolinomialA(M,1,variables[1])*integrand,(variables[1],0,1))
end




#        PolinomialA(M,"t","z")

#This function generates the polinomial A_{t,z}.
#The inputs are a square matrix M and two Char "t".
#    Parameters:
#        M = square Matrix.
#        t,z = characters (a letter inside "").

#    Returns:
#        A dictionary 'components_of_A' whiches keys are the 
#        product of M with the polinomial of z^{n-1}, and the
#        values are elements t^{n}//n.
function PolinomialA(M,t,z)
    s = size(M,1)
    #we need variables to work with
    variable1 = Sym(t)
    variable2 = Sym(z)
    variable1 = conjugate(variable1)
    #generate the components of the polinomial with powers of t
    poli_z = Array{SymPy.Sym}(undef, s)
    for i in 1:s
        poli_z[i] = variable2^(i-1)
    end
    #generate the components of the polinomial with powers of z
    poli_t = Array{SymPy.Sym}(undef, s)
    for i in 1:s
        poli_t[i] = (1 // i)*(variable1^i)
    end
    #multiply matrix times vector
    result = M * poli_z
    #use dot product for the two vectors
    result = dot(poli_t,result)
    #create a dict so that its keys are the components of 'result' and its values are
    #t^n
    #components_of_A = Dict(result[i] => poli_t[i] for i in eachindex(result))
    return result
end




