using BSeries
using SymPy
using LinearAlgebra
using PyCall

# This function generates a polinomial 
#       A_{t,z} = [t,t^2/2,..., t^s/s]*M*[1, z, ..., z^(s-1)]^T
# for a given square matrix M of dimmension s and chars 't' and 'z'.  
function PolinomialA(M,t,z)
    s = size(M,1)
    # we need variables to work with
    variable1 = Sym(t)
    variable2 = Sym(z)
    # conjugate the variable 1, provided that this will be the variable
    # of the left polinomial
    variable1 = conjugate(variable1)
    # generate the components of the polinomial with powers of t
    poli_z = Array{SymPy.Sym}(undef, s)
    for i in 1:s
        poli_z[i] = variable2^(i-1)
    end
    # generate the components of the polinomial with powers of z
    poli_t = Array{SymPy.Sym}(undef, s)
    for i in 1:s
        poli_t[i] = (1 // i)*(variable1^i)
    end
    # multiply matrix times vector
    result = M * poli_z
    #println(result)
    # use dot product for the two vectors
    result = dot(poli_t,result)
    return result
end


"""
    elementary_differentials_csrk(M,tree)

This function calculates the CSRK elementary differential for a given 
square matrix 'M' and a given RootedTree.

"""
function elementary_differentials_csrk(M,rootedtree)
    # we'll work with the level sequence
    tree = rootedtree.level_sequence
    m = maximum(tree)
    l = length(tree)
    # create the variables called 'xi' for 1 <= i <= m
    variables = []
    for i in 1:m
        var_name = "x$i"
        var = Sym(var_name)
        push!(variables, var)
    end
    inverse_counter = l-1
    # stablish initial integrand, which is the rightmost leaf (last node of the level sequence)
    if l > 1
        integrand = integrate(PolinomialA(M,variables[tree[end]-1],variables[tree[end]]),(variables[tree[end]],0,1))
        #println(integrand)
    else
        # if the Rooted Tree is [1] or [], the elementary differential will be 1.
        return 1
    end
    while inverse_counter > 1
        println("flag")
        # println("A_",variables[tree[inverse_counter]-1],variables[tree[inverse_counter]])
        pseudo_integrand = PolinomialA(M,variables[tree[inverse_counter]-1],variables[tree[inverse_counter]])*integrand
        integrand = integrate(pseudo_integrand,(variables[tree[inverse_counter]],0,1))
        # println(integrand)
        inverse_counter -= 1
    end
    #println(integrate(PolinomialA(M,1,variables[1])*integrand,(variables[1],0,1)))
    #multiply for the Basis_Polonimial, i.e. the Polinomial B
    #return the integral with respect to x1.
    return integrate(PolinomialA(M,1,variables[1])*integrand,(variables[1],0,1))
end