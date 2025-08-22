module FerriteHyperelastic
using Ferrite
using Tensors
using LinearAlgebra
using Printf
using GeometryBasics
using Roots

##################### end of packages
export InputStruct
include("InputStruct.jl")
################# end of InputStruc.jl ###############


export assemble_cell_plane_strain!
export assemble_global_plane_strain!
include("utils_plane_strain.jl")
##########################################
##########################################

export solve_lambda3
export assemble_cell_plane_stress!
export assemble_global_plane_stress!
include("utils_plane_stress.jl")
##########################################
##########################################


export assemble_cell_3D!
export assemble_global_3D!
include("utils_threeD.jl")
##########################################
##########################################


export assemble_traction_forces_twoD!
export assemble_traction_forces_threeD!
export run_plane_strain
export run_plane_stress
export run_threeD
export run_fem
include("functions.jl")
##########################################
##########################################
export Faces, Nodes
export to_geometry
export to_boundary
export get_faces
export canonical_quad
export get_boundary_faces
export tet_faces_as_tuples
export canonical_tri
include("utils_plot.jl")


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



end # end of modulus
