module FerriteHyperelastic
using Ferrite
using Tensors
using LinearAlgebra
using Printf
##################### end of packages
export InputStruct
include("InputStruct.jl")
################# end of InputStruc.jl ###############


export elasticity_tensor_plane_stress
export elasticity_tensor_plane_strain
export local_stiffness_linear!
export global_stiffness_linear!


export local_stiffness_nonlinear_plane_strain!
export global_siffness_nonlinear_plane_strain!
export assemble_internal_force_plane_stain!
export assemble_global_internal_force_plane_strain!

export assemble_traction_forces!

export local_stiffness_nonlinear_plane_stress!
export global_siffness_nonlinear_plane_stress!
export assemble_internal_force_plane_stress!
export assemble_global_internal_force_plane_stress!

export run_plane_strain
export run_plane_stress
export run_fem
include("functions.jl")
########## end of functions


end # end of modulus
