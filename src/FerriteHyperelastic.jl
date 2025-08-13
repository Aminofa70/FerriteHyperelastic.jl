module FerriteHyperelastic
using Ferrite
using Tensors
using LinearAlgebra
using Printf
using GeometryBasics

##################### end of packages
export InputStruct
include("InputStruct.jl")
################# end of InputStruc.jl ###############


export elasticity_tensor_plane_stress
export elasticity_tensor_plane_strain
export elasticity_tensor_3D

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


export local_stiffness_nonlinear_3D!
export global_siffness_nonlinear_3D!
export assemble_internal_force_3D!
export assemble_global_internal_force_3D!



export run_plane_strain
export run_plane_stress
export run_3d

export run_fem
include("functions.jl")
########## end of functions

export Faces, Nodes
export to_geometry
export to_boundary
export get_faces
export canonical_quad
export get_boundary_faces
export tet_faces_as_tuples
export canonical_tri
include("utils_plot.jl")

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
