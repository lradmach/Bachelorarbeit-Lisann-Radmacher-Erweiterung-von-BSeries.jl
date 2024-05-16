# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# # Two-Derivative Runge-Kutta Methods

# First, we need to import the packages that will be used in the following code.

import Pkg; Pkg.add("BSeries")
import Pkg; Pkg.add("LinearAlgebra")
import Pkg; Pkg.add("RootedTrees")
# The functions "elementary_weight" and "derivative_weight" also need to be imported to define them new for these methods 
import RootedTrees.elementary_weight
import RootedTrees.derivative_weight

using BSeries
using LinearAlgebra
using RootedTrees

# ## Two-Derivative Runge-Kutta Method as a new AbstractTimeIntegrationMethod
# Note: This is mainly copied from the AdditiveRungeKuttaMethod from Julia package "bseries.jl". 

abstract type AbstractTimeIntegrationMethod end

struct TwoDerivativeRungeKuttaMethod{T, RKm <: AbstractVector{<:RungeKuttaMethod{T}}} <:
       AbstractTimeIntegrationMethod
    rkm::RKm
end

function TwoDerivativeRungeKuttaMethod(rkm) # if not all RK methods use the same eltype
    T = mapreduce(eltype, promote_type, rkm)
    As = map(rk -> T.(rk.A), rkm)
    bs = map(rk -> T.(rk.b), rkm)
    cs = map(rk -> T.(rk.c), rkm)
    TwoDerivativeRungeKuttaMethod(As, bs, cs)
end

function TwoDerivativeRungeKuttaMethod(As, bs, cs = map(A -> vec(sum(A, dims = 2)), As))
    rkm = map(RungeKuttaMethod, As, bs, cs)
    TwoDerivativeRungeKuttaMethod(rkm)
end

Base.eltype(tdrk::TwoDerivativeRungeKuttaMethod{T}) where {T} = T

function Base.show(io::IO, tdrk::TwoDerivativeRungeKuttaMethod)
    print(io, "TwoDerivativeRungeKuttaMethod{", eltype(tdrk), "} with methods\n")
    for (idx, rk) in enumerate(tdrk.rkm)
        print(io, idx, ". ")
        show(io, rk)
    end
end

# ## Elementary Weight of the Two-Derivative Runge-Kutta Method
# As you can see in my bachelor thesis in section ... we need to define the elementary weight completely new for two-derivative Runge-Kutta methods.
# This will be splitted in two parts.

# ### Part One
# Taking a look at section ... in my thesis it is clear that we not only need the 'normal' subtrees of a tree $t$, but also the subtrees of $t$ remaining after removing the tree [1,2] of $t$. Unfortunately, for most of the trees there exist more than one option to remove [1,2] from the tree.
# This brings us to the function "specialsubtrees". It returns a vector of all possibilties of specialsubtrees.

function specialsubtrees(t::RootedTree)
    # first we need an overview about the subtrees of t
    thesubtrees = subtrees(t)
    numberofsubtrees = length(thesubtrees)

    # This will be the vector returning all the possibilities of specialsubtrees
    listofspecialsubtrees = []

    # Using recursion we get the vector of all possibilites of specialsubtees
    for i in 1:numberofsubtrees
        thesubsubtrees = subtrees(thesubtrees[i])
        numberofsubsubtrees = length(thesubsubtrees)

        for j in 1:numberofsubtrees
            if j == i
            elseif j > i
                push!(thesubsubtrees, thesubtrees[j])
            else j < i
                pushfirst!(thesubsubtrees, thesubtrees[i-j])
            end
        
        end

        push!(listofspecialsubtrees, thesubsubtrees)
    end

    return listofspecialsubtrees
end

# ### Part Two
# Now we can start to define the functions "elementary_weight" and "derivative_weight".

function elementary_weight(t::RootedTree, tdrk::TwoDerivativeRungeKuttaMethod)
    b1 = (tdrk.rkm)[1].b
    b2 = (tdrk.rkm)[2].b
    # let a be the elementary weight and n the derivative weight
    # Then the elementary weight is calculated by a(t) = b1 * nu(subtrees(t)) + b2 * nu(specialsubtrees(nu))
    dot(b1, derivative_weight(t, tdrk)) + dot(b2, derivative_weight(t, 0, tdrk))
end

# Since the elementary weight is calculated as written above, we need two different functions for calculating the elementary weight. The first one calculates it using the subtrees of $t$, the second one using the specialsubtrees of $t$.

function derivative_weight(t::RootedTree, tdrk::TwoDerivativeRungeKuttaMethod)
    A1 = (tdrk.rkm)[1].A
    c1 = (tdrk.rkm)[1].c
    A2 = (tdrk.rkm)[2].A
    c2 = (tdrk.rkm)[2].c

    # This vector has the same length like c1 and c2 do but contains only the element 1
    result1 = zero(c1) .+ one(eltype(c1))

    if t == rootedtree(Int64[]) || t == rootedtree([1])
        return zero(c1) .+ one(eltype(c1))

    else
        # Using recursion we calculate the derivative weight as defined in the thesis 
        subtreesoft = subtrees(t)
        numberofsubtreesoft = length(subtreesoft)
        l = 1
        for n in SubtreeIterator(t)
            tmp = A1 * derivative_weight(subtreesoft[l], tdrk) .+ A2 * derivative_weight(subtreesoft[l], 0, tdrk)
            result1 = result1 .* tmp
            l = l + 1
        end
        return result1
    end
end

function derivative_weight(t::RootedTree, a::Int64, tdrk::TwoDerivativeRungeKuttaMethod)
    A1 = (tdrk.rkm)[1].A
    c1 = (tdrk.rkm)[1].c
    A2 = (tdrk.rkm)[2].A
    c2 = (tdrk.rkm)[2].c
    emptyresult = zero(c2)
    
    if t == rootedtree(Int64[])
        return zero(c1) .+ one(eltype(c1))
    else
        
        relevanttreecombinations = specialsubtrees(t)
        number1 = length(relevanttreecombinations)

        # Using recursion we calculate the derivative weight as defined in the thesis
        for k in 1:number1
            relevantsubtreecombinations = relevanttreecombinations[k]
            number2 = length(relevantsubtreecombinations)
            result2 = zero(c1) .+ one(eltype(c1))

            for m in 1:number2
                step = A1 * derivative_weight(relevantsubtreecombinations[m], tdrk) .+ A2 * derivative_weight(relevantsubtreecombinations[m], 0, tdrk)
                result2 = result2 .* step
            end

            emptyresult = emptyresult .+ result2

        end

        return emptyresult
    end
end

# ## B-Series of the Two-Derivative Runge-Kutta Method
# This is mainly copied from the creation of the B-series of RungeKuttaMethod from the Julia package "bseries.jl".

function bseries(tdrk::TwoDerivativeRungeKuttaMethod, order)
    V_tmp = eltype(tdrk)
    if V_tmp <: Integer
        # If people use integer coefficients, they will likely want to have results
        # as exact as possible. However, general terms are not integers. Thus, we
        # use rationals instead.
        V = Rational{V_tmp}
    else
        V = V_tmp
    end
    series = TruncatedBSeries{RootedTree{Int, Vector{Int}}, V}()

    series[rootedtree(Int[])] = one(V)
    for o in 1:order
        for t in RootedTreeIterator(o)
            #only difference here: we are using the elementary_weight for Two-Derivative Runge-Kutta Methods
            series[copy(t)] = elementary_weight(t, tdrk)
        end
    end

    return series
end

function bseries(A::Vector{AbstractMatrix}, b::Vector{AbstractVector}, c::Vector{AbstractVector}, order)
    tdrk = TwoDerivativeRungeKuttaMethod(A, b, c)
    bseries(tdrk, order)
end
