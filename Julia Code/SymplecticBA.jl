# -*- coding: utf-8 -*-
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

# # Checking the symplecticity of a method

# First, we need to import the packages that will be used in the following code.

import Pkg; Pkg.add("BSeries")
import Pkg; Pkg.add("RootedTrees")
import Pkg; Pkg.add("Latexify")

using BSeries
using RootedTrees
using Latexify

# ## The Interesting Part
# Now we can start to define the functions needed to check the symplecticity of a method.
#
# The ideas mainly come from Sundklakk's master thesis and are only changed a little to translate it from python to julia.

# First, we need a function that checks if the condition is satisfied for pair of trees of a specific order.

function SatisfiedForTreesOrder(condition, b, order)
    # Checking the condition for a pair of trees
    for order1 in range(1, order - 1)
        order2 = order - order1
        for tree1 in map(copy, RootedTreeIterator(order1))
            for tree2 in map(copy, RootedTreeIterator(order2))
                if condition(b, tree1, tree2) != true 
                    return false
                end
            end
        end
    end
    return true
end

# Now we define the condition needed for the symplecticity.

function SymplecticityCondition(b::TruncatedBSeries, u, v)
    # Defines the condition needed for symplecticity as mentioned in my thesis
    return b[u∘v] + b[v∘u] == b[u] * b[v]
end

# Finally, we can define the function that returns the symplecticity order of a method.

function Symplecticity(b::TruncatedBSeries, max_order::Int64)
    order = 1

    # Checks if the method is consistent
    if b[rootedtree(Int64[])] != 1
        return 0
        
    else

        # This loops checks the symplecticity condition for a pair of trees of a given order
        for i in 1:max_order
            if in(rootedtree(collect(1:i)), keys(b))
                if SatisfiedForTreesOrder(SymplecticityCondition, b, i) == true                        
                    order = i
                else 
                    break 
                end
            else
                # If the max_order is higher than the order of b the result can be 'wrong' since the coefficients of a(t) of trees
                # with a higher order than order of b are zero
                order = (order, "ATTENTION: max_order higher than order of b. This can affect the symplecticity order. Try again with order of b higher or equal than max_order.")
                break
            end
        end
    end
    return order
end

# # Example
# ## Implicit Midpoint Method
# A good example is the Gauß Method with 2 stages (Implicit Midpoint).
# Since the method is not only pseudo-symplectic but indeed symplectic (means we are having a pseudo-symplecticity order $\infty$, we can choose the max oder as high as we want to since the symplecticity order returned by the function will always be the same as the max oder. 

A = [1//2;;]
b = [1]
c = [1//2]
rk_implicitmidpoint = RungeKuttaMethod(A,b)
bseries_implicitmidpoint = bseries(rk_implicitmidpoint, 3)

Symplecticity(bseries_implicitmidpoint, 6)

# ## Any other explicit method

k = 3
A = [0 0 0;
(8*k-3)//(8k-4) 0 0;
(16k^2-8*k+1)//(2*k*(8k-3)) (2k-1)//(2*k*(8k-3)) 0]
b = [(3k-1)//(8k-3), (-2*(2k-1)^2)//(8k-3), k]
c = [0, (8k-3)//(8k-4), 1]
rkh = RungeKuttaMethod(A,b)
bser = bseries(rkh, 5)


