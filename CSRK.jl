#Packages which would be useful eventually

using BSeries
import SymPy; sp=SymPy; sm = Sym();
using LightGraphs
using Combinatorics
using LinearAlgebra

M = [1 2 0 1 1;
     4 5 6 4 5;
     7 8 9 0 1;
     1 3 4 3 1;
     1 0 2 1 0]

M = [ 1 2 3;
     2 3 4;
     1 0 2]

function Elementary_Weights(M)

end

"""
        PolinomialA(M,t,z)
        
This function generates the polinomial A_{t,z}.
The inputs are a square matrix M and two Char "t".
    Parameters:
        M = square Matrix.
        t,z = characters (a letter inside "").

    Returns:
        A dictionary 'components_of_A' whiches keys are the 
        product of M with the polinomial of z^{n-1}, and the
        values are elements t^{n}//n.
"""
function PolinomialA(M,t,z)
    s = size(M,1)
    #we need variables to work with
    variable1 = Sym(t)
    variable2 = Sym(z)
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
    #result = dot(poli_t,result)
    #create a dict so that its keys are the components of 'result' and its values are
    #t^n
    components_of_A = Dict(result[i] => poli_t[i] for i in eachindex(result))
    return components_of_A
end

"""
This function generates the Polinomial B_zeta for a given 
    square matrix M and a variable zeta = "aletter".
    The code is a simplified version of 'PolinomialA', provided 
    that B is the specific case for which t = 1.
        Parameters:
        M = square Matrix.
        zeta = a char (a letter inside "").

        Returns:
        A dictionary 
"""
function PolinomialB(M,zeta)
    return PolinomialA(M,1,zeta)
end



