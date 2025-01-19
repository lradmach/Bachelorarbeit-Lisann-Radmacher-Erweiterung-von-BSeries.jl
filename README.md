# Bachelorarbeit-Lisann-Radmacher

This Repository is the work that was made during my bachelor thesis. It is an expansion of the "BSeries.jl" package.
There are three topics that were implemented during the last weeks. These topics are symmetry and symplecticity of methods as well as two-derivative Runge-Kutta methods.

The folder "Notebooks" contains three notebooks, one notebook for each topic.
The folder "Julia Code" contains Julia files, one file for each topic.

The implementation for checking the symmetry of a method can be found in the file "SymmetryBA.ipynb" (or "SymmetryBA.jl"). 
In this file you'll find the implementation of the inverse and the adjoint of a given B-Series as well as the function returning up to what order a method is symmetric.

The implementation for checking the symplecticity of a method can be found in the file "SymplecticBA.ipynb" (or "SymplecticBA.jl").
In this file you'll find the implementation of the symplecticity condition, the function which checks the condition for a pair of trees and the function returning up to what order a method is symplectic.

The implementation for two-derivative Runge-Kutta methods can be found in the file "TDRK.ipynb" (or "TDRK.jl").
In this file you'll find the implementation of TDRK methods as a new time integration method, the elementary weight of a TDRK method as well as the B-series of a TDRK method. In addition, there is a bigger example discussed in these files. You can find it in my thesis in the part "Anwendungsbeispiel". The two plots that were made are included in the thesis. You can also see them in the notebook file "TDRKBA.ipynb". 
There is also the option to produce the plots from the shell:
1) Open the "TDRKBA.jl" file
2) include the code
```julia
    include("TDRKBA.jl")
```
3) run the function "plot_firstode()" to produce the plot for the first ODE
```julia
    plot_firstode()
```
4) run the function "plot_secondode()" to produce the plot for the second ODE
```julia
    plot_secondode()
```


The code is procuded with the Julia version 1.9.3.

ABSTRACT

Numerical methods can be used to solve ordinary differential equations. A significant part of these numerical methods are Runge-Kutta methods, which can possess interesting properties. For example, they can be symmetric or symplectic, as extensively discussed in the book Geometric Numerical Integration by Hairer, Lubich, and Wanner (2006). The core of this thesis is formed by B-series, which are based on the concept of rooted trees. These B-series were used throughout this work to analyze the properties of numerical methods using their coefficients, with the master's thesis by H. Sundklakk (2015) serving as the foundation for this analysis.
With the concept of B-series, and especially through the software package BSeries.jl by Ranocha and Ketcheson (2023) for the Julia programming language, along with its extensions developed during this work, an efficient analysis of numerical methods becomes possible.
