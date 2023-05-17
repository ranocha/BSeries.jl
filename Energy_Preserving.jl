# Load the packages we will use.
# These must first be installed using: import Pkg; Pkg.add("package_name")
#import Pkg; Pkg.add("BSeries")
using BSeries
using Latexify  # Only needed for some pretty-printing cells below using `latexify`
import SymPy; sp=SymPy;
using LightGraphs
using Combinatorics
using RootedTrees
using LinearAlgebra

#Parametros iniciales
#stisfies for order <5, and not for 5.
A = [ 0 0 0 
1//3 0 0
-5//48 15//16 0
]
b = [1//10, 1//2, 2//5]


#up to order 3
A = [  
0 0 
2//3 0 
]
b = [1//4, 3//4]

#RK method
A = [0//1  0//1  0//1  0//1
 1//2  0//1  0//1  0//1
 0//1  1//2  0//1  0//1
 0//1  0//1  1//1  0//1]
 b =  [1//6, 1//3, 1//3, 1//6]

#this function checks whether a method is energy Preserving for a given order s
function EnergyPreserving(A,b,s)
    rka = RungeKuttaMethod(A, b)
#generate bseries 
    series_a = modified_equation(bseries(rka, s))
    coefficients = collect(values(series_a))
    atrees = collect(keys(series_a))
# Create an empty vector to store the converted trees
    trees = Vector{Vector{Int}}(undef, length(series_a))
# Convert the trees and store them in the trees vector
    for i in 1:length(series_a)
        #comment: check Rooted tree object
        generate_arrays_from_rootedtree = root_converter(atrees[i])
        if isempty(generate_arrays_from_rootedtree)
            trees[i] = Int[]
        else
            trees[i] = generate_arrays_from_rootedtree
        end
    end
    coefficients = symfact_normalization(coefficients,trees)  
    signal = IsEnergyPreserving(trees,coefficients)  
    if signal == false
        println("Condition Not Satisfied")
    else
        println("Condition Satisfied")
    end
end

function get_t_arrays(a)
    t_dict = Dict{Int, Array}()
    k = a[end]
    for j in 1:k-1
        last_j_occurrence = findlast(x -> x == j, a)
        last_jplus1_occurrence = findlast(x -> x == j+1, a)
        if isnothing(last_j_occurrence) || isnothing(last_jplus1_occurrence)
            t_dict[j] = []
        else
            t_dict[j] = a[last_j_occurrence+1:last_jplus1_occurrence-1]
        end
    end
    return t_dict
end

function modify_t_sub(a)
    t_dict = get_t_arrays(a)
    m = num_leafs(a)
    modified_t_dict = Dict{Int, Vector{Int}}()
    mid_tree = (m+1)/2
    for j in 1:m
        if m % 2 == 1 && j == mid_tree
            modified_t_dict[j] = t_dict[j]
        elseif m % 2 == 1 && j < mid_tree
            new_arr = [n+m-2j+1 for n in t_dict[j]]
            modified_t_dict[j] = new_arr
        elseif m % 2 == 1 && j > mid_tree
            new_arr = [m + n - 2*(j) + 1 for n in t_dict[j]]
            modified_t_dict[j] = new_arr
        elseif m % 2 == 0
            if j <= m/2
                new_arr = [n+m-2j+1 for n in t_dict[j]]
                modified_t_dict[j] = new_arr
            elseif j > m/2
                new_arr = [m + n - 2*(j) + 1  for n in t_dict[j]]
                modified_t_dict[j] = new_arr
            end
        end
    end
    return modified_t_dict
end

function bush_detector(tree)
    bush = false
    if tree != [1]
        l = length(tree)
        bush = true
        if l > 1
            for i in 2:l
                if tree[i] != 2
                    bush = false
                end
            end
        end
    end
    return bush
end
#this function eliminates repeated subarrays
function eliminate_repeated_subarrays(M)
    unique_M = []
    for arr in M
        if !(arr in unique_M)
            push!(unique_M, arr)
        end
    end
    return unique_M
end


#this function swaps the trees
function permuta(a::Vector{Int})
    diccionario = modify_t_sub(a)
    m = num_leafs(a)
    ad_dict = Dict{Int, Vector{Int}}()
    for j in 1:m
        ad_dict[j] = diccionario[m-j+1]
    end
    return ad_dict
end



function adjoint(a::Vector{Int})
    ad_dict = permuta(a)
    m = num_leafs(a)
    adjunto = collect(1:m+1)
    for j in 1:m
        last_j_occurrence = findlast(x -> x == j, adjunto)
        last_jplus1_occurrence = findlast(x -> x == j+1, adjunto)
        adjunto = vcat(adjunto[1:last_j_occurrence], ad_dict[j], adjunto[last_jplus1_occurrence:end])
    end
    adjunto = canonicalarray(adjunto)
    #println("Adjoint Array: ", adjunto)
    return adjunto
end

function root_converter(t::RootedTree{Int64, Vector{Int64}})
    str = string(t) 
    arr_str = match(r"\[(.*?)\]", str).captures[1]
    if isempty(arr_str) || all(iswhitespace, arr_str)
        return Int[]
    end
    arr = parse.(Int, split(arr_str, ","))
    return arr
end
#this function reorders the tree in the same way as the Bseries output does
function canonicalarray(tree)
    t = rootedtree(tree)
    trees = root_converter(t)
    return trees
end

#this functions receives
function indexator(trees,arbol)
    variable = 0
    for i in 1:length(trees)
        if trees[i] == arbol
            variable = i
        end
    end
    return variable
end

iswhitespace(c::Char) = isspace(c)

function BMinus(array)
    m = 1
    for i in array
        if i == 2
            m = m + 1
        end
    end
    contador = 0
    #aqui almacenamos los 
    auxiliar = Array{Int64}(undef, m)
    longitud = length(array)
    auxiliar[m] = longitud
    #numero de subtrees
    n = m-1
    subtrees = Vector{Any}(undef, n)
    for i in 2:longitud
        if array[i] == 2
            contador += 1
            auxiliar[contador] = i
            if contador == n
                break
            end
        end
    end
    for i in 1:n
        subtrees[i] = array[auxiliar[i]:(auxiliar[i+1]-1)]
    end
    push!(subtrees[n],array[longitud])
    return subtrees
end

#println(BMinus(array))


#this function returns a Hamiltonian H matriz of nxn


              
#funcion generadora de arboles simetricos
function symmetries(tree)
    subtrees = BMinus(tree)
    #first order symmetries
    perms = collect(permutations(subtrees))
    l = length(perms)
    for i in 1:l
        perms[i] = reduce(vcat, perms[i])
    end
    result = add_one_to_left(perms)
    return(result)
end

#generates all equivalent trees
function equivalent_trees(array)
    tree = rootedtree(array)
    l = length(array)
    superarray = get_permutations(array)
    lperm = length(superarray)
    for i in reverse(1:lperm) 
        if superarray[i][1] != 1
            splice!(superarray, i)
        end
    end
    lperm = length(superarray)
    flag = false
    for i in reverse(1:lperm)
        for j in 1:(l-1)
            if (superarray[i][j+1] - superarray[i][j]) > 1
                flag = true
            end
        end
        if flag == true
            splice!(superarray, i)
        end
        flag = false
    end
    lperm = length(superarray)
    for i in reverse(1:lperm)
        arbol = rootedtree(superarray[i])
        bandera = arbol == tree
        if bandera == false
            splice!(superarray, i)
        end
    end
    bandera = false
    lperm = length(superarray)
    for i in reverse(1:lperm)
        for j in 1:(i-1)
            if superarray[i] == superarray[j]
                bandera = true
            end
        end 
        if bandera == true
            splice!(superarray, i)
        end
        bandera = false
    end

    return superarray
end

function get_permutations(arr)
    perms = collect(permutations(arr))
    return perms
end

function num_leafs(a)
    k = a[end]
    return k - 1
end

function add_one_to_left(arrays)
    for i in 1:length(arrays)
        arrays[i] = vcat([1], arrays[i])
    end
    return arrays
end

#this function multplies the coefficient for its symmetry factor
function symfact_normalization(coef,thetrees)
    l = length(coef)
    for i in 1:l
        factor = 0
        factor = RootedTrees.symmetry(RootedTree(thetrees[i]))
        coef[i] = coef[i]*(1//factor)
    end
    return coef
end

#this function tells up to what order a method is Energy Preserving
function OrderMethod(A,b)
    s = 1
    flag = false
    rka = RungeKuttaMethod(A, b)
    while flag == false 
#generate bseries 
        series_a = modified_equation(bseries(rka, s))
        coefficients = collect(values(series_a))
        atrees = collect(keys(series_a))
# Create an empty vector to store the converted trees
        trees = Vector{Vector{Int}}(undef, length(series_a))
# Convert the trees and store them in the trees vector
        for i in 1:length(series_a)
            array_from_RTree = root_converter(atrees[i])
            if isempty(array_from_RTree)
                trees[i] = Int[]
            else
                trees[i] = array_from_RTree
            end
        end
        coefficients = symfact_normalization(coefficients,trees)  
        #println(trees)
        #println(coefficients)
        if IsEnergyPreserving(trees,coefficients) == true
            flag = true
        end
        s = s + 1
    end
    println("Energy Preserving for order < ", s-1)
end

#this function checks if an AVF method is energy Preserving for a given order s
function EnergyPreservingAVF(s) 
    series = bseries(s) do t, series
        if order(t) in (0, 1)
            return 1 // 1
        else
            v = 1 // 1
            n = 0
            for subtree in SubtreeIterator(t)
                v *= series[subtree]
                n += 1
            end
            return v / (n + 1)
        end
    end
#generate bseries 
    series_a = modified_equation(series)
    coefficients = collect(values(series_a))
    atrees = collect(keys(series_a))
    trees = Vector{Vector{Int}}(undef, length(series_a))
# Convert the trees and store them in the trees vector
    for i in 1:length(series_a)
        array_from_RTree = root_converter(atrees[i])
       if isempty(array_from_RTree)
            trees[i] = Int[]
        else
            trees[i] = array_from_RTree
        end
    end
    coefficients = symfact_normalization(coefficients,trees)  
    signal = IsEnergyPreserving(trees,coefficients)  
    if signal == false
        println("Condition Not Satisfied")
    else
        println("Condition Satisfied")
    end
    #series_a
end

#do not run this function 
#runs infty
function OrderAVF() 
    s = 1
    flag = false
    while flag == false 
        series = bseries(s) do t, series
            if order(t) in (0, 1)
                return 1 // 1
            else
                v = 1 // 1
                n = 0
                for subtree in SubtreeIterator(t)
                    v *= series[subtree]
                    n += 1
                end
                return v / (n + 1)
            end
        end
#generate bseries 
        series_a = modified_equation(series)
        
        coefficients = collect(values(series_a))
        atrees = collect(keys(series_a))
# Create an empty vector to store the converted trees
        trees = Vector{Vector{Int}}(undef, length(series_a))
# Convert the trees and store them in the trees vector
        for i in 1:length(series_a)
            array_from_RTree = root_converter(atrees[i])
            if isempty(array_from_RTree)
                trees[i] = Int[]
            else
                trees[i] = array_from_RTree
            end
        end
        coefficients = symfact_normalization(coefficients,trees)  
        #println(trees)
        #println(coefficients)
        if IsEnergyPreserving(trees,coefficients) == false
            flag = true
        end
        s = s + 1
    end
    #println("bandera", s)
    println("Energy Preserving for order < ", s-1)    
end


 
function IsEnergyPreserving(trees, coefficients)
    #obtain the adjoint and check if it exists
    numero_de_coef = length(trees)
    MatrizEP = Vector{Int64}[]
    signal = false
    #bigflag = false
    #bandera_array = false
    for t in trees
        if !isempty(t) 
            indice_mayor = indexator(trees,t)
            #println(indice_mayor)
            bushflag = bush_detector(t)
            if bushflag == true
                if coefficients[indice_mayor] != 0
                    signal = true
                end
                #println("bush")
            else
                equiv_set = equivalent_trees(t)
                bandera_ad = false
                for arbol in equiv_set
                #j-th canonical vector
                    ej = zeros(Int64, numero_de_coef)
                    ej[indice_mayor] = 1
                #continue
                    m = num_leafs(arbol)
                    t_ad = adjoint(arbol)
                #revisamos si hay self-adjoint
                    if (rootedtree(t_ad) == rootedtree(arbol))
                        bandera_ad = true
                    end
                    ek = zeros(Int64, numero_de_coef) 
                    if t_ad in trees 
                    #bandera_array = true
                    #k-th canonical vector
                        ek[indexator(trees,t_ad)] = 1*(-1)^m
                    #println(ek)
                        ej = ej + ek
                        #println(ej)
                        push!(MatrizEP, ej)
                    end
                end
            end
        end
    end
    #eliminamos las columnas vacías
    filter!(x -> any(y -> y != 0, x), MatrizEP)
    #eliminamos las columnas repetidas
    MatrizEP = eliminate_repeated_subarrays(MatrizEP)
    #println(MatrizEP)
    #println("Longitud de M: ", length(MatrizEP))
    M = hcat(MatrizEP...)
    rank_M = rank(M)
    rank_MV = rank([M coefficients])
    #println(rank_M == rank_MV)
    #X = find_vector_X(M,coefficients)
    #println(X)
    println("All the trees have been checked:")
    if signal == true
        println("Condition Not Satisfied")
    else
        println("Condition Satisfied")
    end
    return rank_M == rank_MV
end
