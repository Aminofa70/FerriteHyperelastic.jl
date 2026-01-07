module FerriteHyperelastic
############################
using Ferrite
using Tensors
using LinearAlgebra
using Printf
using GeometryBasics
using Roots
using IterativeSolvers
using BlockArrays
using SparseArrays
using Statistics
##################### end of packages



##########################################
##########################################

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
##########################################
##########################################

##########################################
##########################################
end # end of modulus
