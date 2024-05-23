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

# # Checking the symmetry of methods

# First, I need to import the packages that will be used in the following code.

import Pkg; Pkg.add("BSeries")
import Pkg; Pkg.add("RootedTrees")

using BSeries
using RootedTrees

# ## Antipode
# Since we know that the inverse element in the Butcher group is defined as $a^{-1}=a \circ S$, where S is the so-called antipode, we need to define this antipode first.

function antipode(t::RootedTree)
    # The all_splittings function from "RootedTrees.jl" can be used for the coproduct
    coproduct = all_splittings(t)

    # Collecting the trees that are removed from the tree
    forests = coproduct[1]

    # Collecting the subtree that remains after removing the trees
    subtrees = coproduct[2]

    # Number of tensor products 
    number = length(forests)

    # The function antipode returns a linear combination in a vector. This vector contains two vectors:
    # The first vector (coeffs) returns the coefficients
    coeffs = [-1]
    # The second vector (trees) returns the trees
    trees = [[t]]

    # only going from 2:(number-1) since the two tensor products containing the empty tree are not relevant here
    for i in 2:(number-1)

        leftovertree = subtrees[i]

        # This loop pushes the trees in the "trees" vector
        for m in antipode(leftovertree)[2]
            empty2 = []
            #length = length(m)

            for tree in forests[i]
            push!(empty2, tree)
            end
            
            for n in m
                push!(empty2, n)
            end
            
            push!(trees, empty2)
        end

        # This loop pushes the coefficients to the "coeffs" vector
        for j in antipode(leftovertree)[1]
            push!(coeffs, -1*j)
        end
      
    end   
        
    # Returns the linear combination as a vector: that means we sum over all coeffs[i]*trees[i]
    return [coeffs, trees]
end

# Now I am using the antipode for the inverse.

# +
function inverse(b::TruncatedBSeries, t::RootedTree)
    result = 0

    # First calculating the antipode of the given tree t
    antipode_tree = antipode(t)

    # Now calculating the inverse: a^(-1) = a * S
    k = 1
    for i in antipode_tree[1]
        btw = i
        for j in antipode_tree[2][k]
            btw *= b[j]
        end
        k += 1
        result += btw
    end
    return result
    
end
# -

# Finally, I can define the adjoint of a given B-series. The function "adjoint" returns the time-reversed method of the inverse B-series.

function adjoint(b::TruncatedBSeries)
    # Returns a new B-Series: The Adjoint B-Series
    series = TruncatedBSeries{RootedTree{Int, Vector{Int}}, Rational{Int64}}()

    series[rootedtree(Int[])] = one(Rational{Int64})
    orderofb = order(b)
    for o in 1:orderofb
        for t in RootedTreeIterator(o)
            # Changes the sign of a(t) whenever |t| is odd
            series[copy(t)] = (-1)^(order(t)) * inverse(b,t)
        end
    end
    return series
end

# ## Checking up to what order a method if symmetric

function SymmetricOrder(b::TruncatedBSeries)
    # Now checking the order to which the method is symmetric. This works similar to "order_of_accuracy" from "BSeries.jl"
    adj = adjoint(b)

    for o in 0:order(b)
        # Iterate over all rooted trees used as keys in `series`
        # of a given order `o`.
        for t in RootedTreeIterator(o)
            if (b[t] - adj[t]) != 0
                return order(t) - 1
            end
        end
    end

    return order(b)
end

# ## Checking up to what order a method is adjoint to another method

function AdjointOrder(b::TruncatedBSeries, c::TruncatedBSeries)
    adj = adjoint(b)

    # Now checking the order to which the method is adjoint to another method. This works similar to "SymmetricOrder".
    for o in 0:order(b)
        # Iterate over all rooted trees used as keys in `series`
        # of a given order `o`.
        for t in RootedTreeIterator(o)
            if (c[t] - adj[t]) != 0
                return order(t) - 1
            end
        end
    end

    return order(b)
end

# ## Examples for symmetric methods

# ### Implicit Midpoint
# The implicit midpoint method is symmetric. Thus, the symmetric order will be always the same as the order we choose for the B-series.

A_im = [1//2 ;;]
b_im = [1]
c_im = [1//2]
rk_im = RungeKuttaMethod(A_im, b_im, c_im)
bseries_im = bseries(rk_im, 6);

SymmetricOrder(bseries_im)

# ### Explicit Midpoint
# The Explicit Midpoint is only symmetric up to order 3.

A_em = [0 0; 1//2 0]
b_em = [0, 1]
c_em = [0, 1//2]
rk_em = RungeKuttaMethod(A_em, b_em, c_em)
bseries_em = bseries(rk_em, 6);

order_of_accuracy(bseries_em)

SymmetricOrder(bseries_em)

# ## Examples for methods that are adjoint to another method

# ### Implicit and Explicit Euler

A_ex_euler = [0 ;;]
b_ex_euler = [1]
c_ex_euler = [0]
rk_ex_euler = RungeKuttaMethod(A_ex_euler, b_ex_euler, b_ex_euler)
bseries_ex_euler = bseries(rk_ex_euler, 5);

A_im_euler = [1 ;;]
b_im_euler = [1]
c_im_euler = [1]
rk_im_euler = RungeKuttaMethod(A_im_euler, b_im_euler, b_im_euler)
bseries_im_euler = bseries(rk_im_euler, 5);

# The implicit Euler method is the adjoint of the explicit Euler method.

AdjointOrder(bseries_im_euler, bseries_ex_euler)
