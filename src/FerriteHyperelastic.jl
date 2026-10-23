module FerriteHyperelastic
## the packages calling
using Ferrite , Tensors, LinearAlgebra, Printf
using GeometryBasics , Roots, IterativeSolvers , BlockArrays
using SparseArrays, Statistics, GLMakie, CairoMakie
using ComodoFerrite , FESolvers,  HyperData , ForwardDiff
using NonlinearSolve
## export the package to be used in the FerriteHyperelastic 
export Ferrite, Tensors, Printf 
export GeometryBasics, Roots, IterativeSolvers
export BlockArrays, SparseArrays, Statistics
export LinearAlgebra, GLMakie, CairoMakie
export ComodoFerrite, FEBio, FESolvers 
export HyperData, ForwardDiff, NonlinearSolve
## export the functions for the curve-fitting
export calc_init_shear_mod
export calc_mat_constants_hyper_biaxial
export calc_mat_constants_hyper_ogden_linear
export calc_residual_and_jacobian
export calc_mat_constants_hyper_ogden
export calc_mat_constants_hyper_shear
export calc_mat_constants_hyper_uniaxial
export drucker_stability_biaxial
export drucker_stability_hyper
export drucker_stability_ogden
export drucker_stability_shear
export drucker_stability_uniaxial
export solver_constants_hyper
include("utils_fitting.jl")

end # end of modulus
