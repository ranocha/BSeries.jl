#This code is based on the Theorem 2 of the paper "Energy-Preserving Integrators and the Structure of B-series" (link: https://link.springer.com/article/10.1007/s10208-010-9073-1).
#Functions to run: EnergyPreserving(A,b,s), EnergyPreservingAVF(s), BSeries_Energy_Preserving(series)
# Load the packages we will use.
using BSeries
import SymPy; sp=SymPy;
using Combinatorics
#using RootedTrees
using LinearAlgebra
import RootedTrees; rt=RootedTrees


#this function checks whether a method is energy Preserving for a given order s
function EnergyPreserving(A,b,s)
    rka = RungeKuttaMethod(A, b)
#generate bseries 
    series_a = modified_equation(bseries(rka, s))
    #save all the coefficients in an array
    coefficients = collect(values(series_a))
    #save all the RootedTrees in another array: 
    #we need only the level sequence
    atrees = collect(keys(series_a))
# Create an empty vector to store the converted trees into arrays
    trees = Vector{Vector{Int}}(undef, length(series_a))
# Convert the trees and store them in the 'trees' vector
    for i in 1:length(series_a)
        levelsequence = atrees[i].level_sequence
        if isempty(levelsequence)
            trees[i] = Int[]
        else
            trees[i] = levelsequence
        end
    end
    #normalize the coefficients multiplying by the symmetry factor 
    coefficients = symfact_normalization(coefficients,trees) 
    #check if it is energy Preserving 
    signal = IsEnergyPreserving(trees,coefficients)  
    return signal
end

#this function checks if the modified_equation of a BSeries is Energy Preserving or not
function is_energy_preserving(series)
    series_a = modified_equation(series)
    #save all the coefficients in an array
    coefficients = collect(values(series_a))
    #save all the RootedTrees in another array: 
    #we need only the level sequence
    atrees = collect(keys(series_a))
# Create an empty vector to store the converted trees into arrays
    trees = Vector{Vector{Int}}(undef, length(series_a))
# Convert the trees and store them in the 'trees' vector
    for i in 1:length(series_a)
        levelsequence = atrees[i].level_sequence
        if isempty(levelsequence)
            trees[i] = Int[]
        else
            trees[i] = levelsequence
        end
    end
    #normalize the coefficients multiplying by the symmetry factor 
    coefficients = symfact_normalization(coefficients,trees) 
    #check if it is energy Preserving 
    signal = IsEnergyPreserving(trees,coefficients)  
    return signal
end

function get_leafs(a)
    t_dict = Dict{Int, Array}()
    #we need to save the final nnumber in the level_sequence because this is the final leaf of the spine
    k = a[end]
    #we need to look for the last last_j_occurrence of every integer in [1,k-1]
    for j in 1:k-1
        last_j_occurrence = findlast(x -> x == j, a)
        last_jplus1_occurrence = findlast(x -> x == j+1, a)
        #consider the empty leafs
        if isnothing(last_j_occurrence) || isnothing(last_jplus1_occurrence)
            t_dict[j] = []
        else
            t_dict[j] = a[last_j_occurrence+1:last_jplus1_occurrence-1]
        end
    end
    return t_dict
end

#it is not enough to generate the leafs and swap them. The level_sequences must be modified 
#with corrections to the numbers inside: the numbers will decrease if the leaf is moved to a 
#lower position, and they will increase if they move to an upper position.
function modify_t_sub(a)
    #we obtain the leafs via 'get_leafs'
    #save them in 't_dict'
    t_dict = get_leafs(a)
    m = num_leafs(a)
    #create another dict for the modified indexes
    modified_t_dict = Dict{Int, Vector{Int}}()
    mid_tree = (m+1)/2
    #we check if the number of leafs is odd:
    #in that case, the middle one remains the same
    for j in 1:m
        if m % 2 == 1 && j == mid_tree
            modified_t_dict[j] = t_dict[j]
            #now, go for the odd m case:
            #if the original leaf is low (with respect to the middle position), we use the formula
            # n+m-2j+1 for every number in the level_sequence
        elseif m % 2 == 1 && j < mid_tree
            new_arr = [n+m-2j+1 for n in t_dict[j]]
            modified_t_dict[j] = new_arr
            #if it is high, use m + n - 2*(j) + 1
        elseif m % 2 == 1 && j > mid_tree
            new_arr = [m + n - 2*(j) + 1 for n in t_dict[j]]
            modified_t_dict[j] = new_arr
            #even m case
        elseif m % 2 == 0
            #for low: n+m-2j+1
            if j <= m/2
                new_arr = [n+m-2j+1 for n in t_dict[j]]
                modified_t_dict[j] = new_arr
                #for high: m + n - 2*(j) + 1
            elseif j > m/2
                new_arr = [m + n - 2*(j) + 1  for n in t_dict[j]]
                modified_t_dict[j] = new_arr
            end
        end
    end
    return modified_t_dict
end

#this functions checks if the level_sequence is a bush
function bush_detector(tree)
    bush = false
    if tree != [1]
        l = length(tree)
        bush = true
        initial_index = tree[1]
        initial_plus_one = initial_index + 1
        if l > 1
            for i in 2:l
                if tree[i] != initial_plus_one
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
        if arr in unique_M
        else
            push!(unique_M, arr)
        end
    end
    return unique_M
end

#this function swaps the trees
function permuta(a::Vector{Int})
    #generate all the leafs via 'modify_t_sub'
    diccionario = modify_t_sub(a)
    m = num_leafs(a)
    ad_dict = Dict{Int, Vector{Int}}()
    #we swap every leaf j for the (m-j+1)
    for j in 1:m
        ad_dict[j] = diccionario[m-j+1]
    end
    return ad_dict
end


#this function returns the rightmost_energy_preserving_tree level sequence for a given tree (the input also in level-sequence form)
#the adoint is calculated with respect to the right-most spine
function rightmost_energy_preserving_tree(a::Vector{Int})
    #we want to generate all the leafs with respect to the rightmost spine
    ad_dict = permuta(a)
    #we obtain the number of leafs the right-most spine has
    #the Theorem 2 in the article requires to know if m is odd or even
    m = num_leafs(a)
    #we create an array from 1 to m plus another node for the final leaf of the rightmost spine
    adjunto = collect(1:m+1)
    #then, we insert every level_sequence from ad_dict
    for j in 1:m
        last_j_occurrence = findlast(x -> x == j, adjunto)
        last_jplus1_occurrence = findlast(x -> x == j+1, adjunto)
        adjunto = vcat(adjunto[1:last_j_occurrence], ad_dict[j], adjunto[last_jplus1_occurrence:end])
    end
    adjunto = canonicalarray(adjunto)
    #println("Adjoint Array: ", adjunto)
    return adjunto
end

#this function rearranges the level sequence in the same way as the Bseries output does
#
function canonicalarray(tree)
    t = rootedtree(tree)
    trees = t.level_sequence
    return trees
end

#this functions receives as inputs an array-of-arrays and one of its elements.
#the output is the index of this array 'onetree' 
function indexator(trees,onetree)
    theindex = 0
    l = length(trees)
    for i in 1:l
        if trees[i] == onetree
            theindex = i
        end
    end
    return theindex
end

iswhitespace(c::Char) = isspace(c)
#this function is the BMinus operator as described in B-Series papers
function BMinus(array)
    m = 1
    for i in array
        if i == 2
            m = m + 1
        end
    end
    counter = 0
    #we save the trees here
    auxiliar = Array{Int64}(undef, m)
    l = length(array)
    auxiliar[m] = l
    #the number of subtrees
    n = m-1
    subtrees = Vector{Any}(undef, n)
    for i in 2:l
        if array[i] == 2
            counter += 1
            auxiliar[counter] = i
            if counter == n
                break
            end
        end
    end
    for i in 1:n
        subtrees[i] = array[auxiliar[i]:(auxiliar[i+1]-1)]
    end
    push!(subtrees[n],array[l])
    return subtrees
end
         

#generates all equivalent trees
function equivalent_trees(array)
    tree = rootedtree(array)
    l = length(array)
    #we get the permutations of the array
    superarray = get_permutations(array)
    lperm = length(superarray)
    #eliminate all the arrays except for those for which the first element is equal to 1
    for i in reverse(1:lperm) 
        if superarray[i][1] != 1
            splice!(superarray, i)
        end
    end
    lperm = length(superarray)
    #now, check if the levelsequence is the level sequence represents a RootedTree:
    #the condition is that each number m cannot be followed by a number n such that
    #n-m > 1. If that's the case, eliminate this array with 'splice!'
    RTree_array = true
    for i in reverse(1:lperm)
        for j in 1:(l-1)
            if (superarray[i][j+1] - superarray[i][j]) > 1
                RTree_array = false
            end
        end
        if RTree_array == false
            splice!(superarray, i)
        end
        RTree_array = true
    end
    lperm = length(superarray)
    for i in reverse(1:lperm)
        onetree = rootedtree(superarray[i])
        isthesame_tree = onetree == tree
        if isthesame_tree == false
            splice!(superarray, i)
        end
    end
    repeated_tree = false
    lperm = length(superarray)
    #we delete any repeated tree
    for i in reverse(1:lperm)
        for j in 1:(i-1)
            if superarray[i] == superarray[j]
                repeated_tree = true
            end
        end 
        if repeated_tree == true
            splice!(superarray, i)
        end
        repeated_tree = false
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
    l = length(arrays)
    for i in 1:l
        arrays[i] = vcat([1], arrays[i])
    end
    return arrays
end

#this function multplies the coefficient for its symmetry factor
function symfact_normalization(coef,thetrees)
    l = length(coef)
    for i in 1:l
        factor = 0
        #because of the librarys we are using, some packages are repeated
        #Then, we specify that the symmetry function comes from RootedTrees
        factor = rt.symmetry(RootedTree(thetrees[i]))
        coef[i] = coef[i]*(1//factor)
    end
    return coef
end

#this function tells up to what order a method is Energy Preserving
function OrderMethod(A,b)
    s = 1
    energy_preserving = false
    rka = RungeKuttaMethod(A, b)
    while energy_preserving == false 
#generate bseries 
        series_a = modified_equation(bseries(rka, s))
        coefficients = collect(values(series_a))
        atrees = collect(keys(series_a))
# Create an empty vector to store the converted trees
        trees = Vector{Vector{Int}}(undef, length(series_a))
# Convert the trees and store them in the trees vector
        for i in 1:length(series_a)
            array_from_RTree = atrees[i].level_sequence
            if isempty(array_from_RTree)
                trees[i] = Int[]
            else
                trees[i] = array_from_RTree
            end
        end
        coefficients = symfact_normalization(coefficients,trees)  
        if IsEnergyPreserving(trees,coefficients) == false
            energy_preserving = true
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
        array_from_RTree = atrees[i].level_sequence
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
    condition_satisfied = false
    while condition_satisfied == false 
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
            array_from_RTree = atrees[i].level_sequence
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
            condition_satisfied = true
        end
        s = s + 1
    end
    #println("bandera", s)
    println("Energy Preserving for order < ", s-1)    
end

function IsEnergyPreserving(trees, coefficients)
    #for every tree, obtain the adjoint and check if it exists
    length_coeff = length(trees)
    #provided that a tree can have several adjoints, for every tree t we need to check if there is a linear combination of the coefficients 
    #of its adjoints such that this linear combination is equal to the coefficient of t. 
    #For this purpose, we cretae a MatrizEP to store all the information
    MatrizEP = Vector{Int64}[]
    error_searcher = false
    for t in trees
        if !isempty(t) 
            #save the index of the level_sequence: this will be used for creating the matrix EP
            highindex = indexator(trees,t)
            #check if the level_sequence corresponds to a bush
            bushflag = bush_detector(t)
            #if it does, then the coefficient must be equal to zero.
            if bushflag == true
                if coefficients[highindex] != 0
                    error_searcher = true
                end
            else
                #if the tree is not a bush, then generate all the equivalent trees
                equiv_set = equivalent_trees(t)
                #this flag checks if the tree is self-adjoint
                #flag_ad = false
                for onetree in equiv_set
                #j-th canonical vector
                    ej = zeros(Int64, length_coeff)
                    ej[highindex] = 1
                #continue
                    m = num_leafs(onetree)
                    t_ad = rightmost_energy_preserving_tree(onetree)
                #check if the tree t is self-adjoint
                    #if (rootedtree(t_ad) == rootedtree(onetree))
                    #    flag_ad = true
                    #end
                    ek = zeros(Int64, length_coeff) 
                    #check if an rightmost_energy_preserving_tree is in the set of trees
                    if t_ad in trees 
                    #generate a k-th canonical vector
                        ek[indexator(trees,t_ad)] = 1*(-1)^m
                    #println(ek)
                    #sum ej + ek and push it into MatrizEP
                        ej = ej + ek
                        #println(ej)
                        push!(MatrizEP, ej)
                    end
                end
            end
        end
    end
    #we filter the empty columns (those full of zeros)
    filter!(x -> any(y -> y != 0, x), MatrizEP)
    #we filter repeated columns
    MatrizEP = eliminate_repeated_subarrays(MatrizEP)
    #println(MatrizEP)
    #println("Longitud de M: ", length(MatrizEP))
    #because the components of MatrizEP is supposed to be columns, we traspose the matrix
    M = hcat(MatrizEP...)
    rank_M = rank(M)
    #we also create an extended matrix for which we append the vector of coefficients
    rank_MV = rank([M coefficients])
    #println(rank_M == rank_MV)
    #X = find_vector_X(M,coefficients)
    #println(X)
    #println("All the trees have been checked:")
    if error_searcher == true
        result = false
    else
        result = rank_M == rank_MV
    end
    #if the rank of M is equal to the rank of the extended MV, then the system is energy-Preserving
    return result
end
