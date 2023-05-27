#Packages which would be useful eventually

using BSeries
import SymPy; sp=SymPy;
using LightGraphs
using Combinatorics
using LinearAlgebra

M = [1 2 0 1 1;
     4 5 6 4 5;
     7 8 9 0 1;
     1 3 4 3 1;
     1 0 2 1 0]

function Elementary_Weights(M)
    Polinomial(M,t,z)

    tau_array = 
    for i in 1:s

    end
    
end

#this function generates the polinomial A_{t,z}
function PolinomialAtz(M,t,z)
    s = size(M,1)
    #we need variables to work with
    variable1 = Sym("t")
    variable2 = Sym("z")
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
    return result
end



