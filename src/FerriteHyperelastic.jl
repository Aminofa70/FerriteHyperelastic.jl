module FerriteHyperelastic
using Ferrite
using Tensors
##################### end of packages
export InputStruc
include("InputStruc.jl")
################# end of InputStruc.jl ###############


export elasticity_tensor_plane_stress
export elasticity_tensor_plane_strain
export local_stiffness_linear!
export global_stiffness_linear!
export local_stiffness_nonlinear!
export global_siffness_nonlinear!
export assemble_internal_force!
export assemble_global_internal_force!
export assemble_traction_forces!
export run_fem
include("functions.jl")
########## end of functions


end # end of modulus
