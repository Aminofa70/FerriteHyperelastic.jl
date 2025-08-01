"""
elasticity_tensor_plane_stress(E, ν)
Elasticity tensor for plane stress deformation
"""
function elasticity_tensor_plane_stress(E, ν)
    C_voigt = E / (1 - ν^2) * [
        1.0 ν 0.0;
        ν 1.0 0.0;
        0.0 0.0 (1-ν)/2
    ]
    return fromvoigt(SymmetricTensor{4,2}, C_voigt)
end
######################## end of elasticity_tensor_plane_stress ########################
"""
elasticity_tensor_plane_strain(E, ν)

Elasticity tensor for plane strain deformation
"""
function elasticity_tensor_plane_strain(E, ν)
    C_voigt = E / ((1 + ν) * (1 - 2ν)) * [
        1-ν ν 0.0;
        ν 1-ν 0.0;
        0.0 0.0 (1-2ν)/2
    ]
    return fromvoigt(SymmetricTensor{4,2}, C_voigt)
end
########################## end of elasticity_tensor_plane_strain #################
"""
Elasticity tensor for 3D deformation
"""
function elasticity_tensor_3D(E, ν)
    C_voigt = E / ((1 + ν) * (1 - 2 * ν)) * [
        1-ν ν ν 0 0 0;
        ν 1-ν ν 0 0 0;
        ν ν 1-ν 0 0 0;
        0 0 0 (1-2*ν)/2 0 0;
        0 0 0 0 (1-2*ν)/2 0;
        0 0 0 0 0 (1-2*ν)/2
    ]
    return fromvoigt(SymmetricTensor{4,3}, C_voigt)
end
########################## end of elasticity_tensor_3D #################
"""
local_stiffness_linear!(ke, cellvalues, C )

Local Stiffness Matrix
Create local stiffness for linear elastic part
"""
function local_stiffness_linear!(ke, cellvalues, C)
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:getnbasefunctions(cellvalues)
            ∇Nᵢ = shape_gradient(cellvalues, q_point, i)
            for j in 1:getnbasefunctions(cellvalues)
                ∇ˢʸᵐNⱼ = shape_symmetric_gradient(cellvalues, q_point, j)
                ke[i, j] += (∇Nᵢ ⊡ C ⊡ ∇ˢʸᵐNⱼ) * dΩ
            end
        end
    end
    return ke
end
########################### end of local_stiffness_linear!#################
"""
global_stiffness_linear!(K_linear, dh, cellvalues, C)

Global Stiffness Matrix->Assembled from all local stiffness matrices.
"""
function global_stiffness_linear!(K_linear, dh, cellvalues, C)
    n_basefuncs = getnbasefunctions(cellvalues)
    ke = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K_linear)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(ke, 0.0)
        local_stiffness_linear!(ke, cellvalues, C)
        assemble!(assembler, celldofs(cell), ke)
    end
    return K_linear
end
######################## end of global_stiffness_linear! ##################

"""
local_stiffness_nonlinear_plane_strain!(ke_n, cell, cv, input::InputStruct, ue_global)
Local Stiffness Matrix for nonlinear part: plane strain
"""
function local_stiffness_nonlinear_plane_strain!(ke_n, cell, cv, input::InputStruct, ue_global)
    reinit!(cv, cell)
    fill!(ke_n, 0.0)
    ndofs = getnbasefunctions(cv)
    cell_dofs = celldofs(cell)
    ue_local = view(ue_global, cell_dofs)
    material = input.material
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇u2d = function_gradient(cv, qp, ue_local)
        F = Tensor{2,3,Float64}((
            1.0 + ∇u2d[1, 1], ∇u2d[1, 2], 0.0,
            ∇u2d[2, 1], 1.0 + ∇u2d[2, 2], 0.0,
            0.0, 0.0, 1.0
        ))
        C = tdot(F)
        S, ∂S∂C = material(C)
        I = one(S)
        ∂P∂F = otimesu(I, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(F', I)
        for i in 1:ndofs
            ∇δui2d = shape_gradient(cv, qp, i)
            ∇δui = Tensor{2,3,Float64}((
                ∇δui2d[1, 1], ∇δui2d[1, 2], 0.0,
                ∇δui2d[2, 1], ∇δui2d[2, 2], 0.0,
                0.0, 0.0, 0.0
            ))
            ∇δui_∂P∂F = ∇δui ⊡ ∂P∂F
            for j in 1:ndofs
                ∇δuj2d = shape_gradient(cv, qp, j)
                ∇δuj = Tensor{2,3,Float64}((
                    ∇δuj2d[1, 1], ∇δuj2d[1, 2], 0.0,
                    ∇δuj2d[2, 1], ∇δuj2d[2, 2], 0.0,
                    0.0, 0.0, 0.0
                ))
                ke_n[i, j] += (∇δui_∂P∂F ⊡ ∇δuj) * dΩ
            end
        end
    end
    return
end
############################ end of local_stiffness_nonlinear! ###########################
"""
global_siffness_nonlinear_plane_strain!(K_nonlinear, dh, cellvalues, input::InputStruct, ue)

"""
function global_siffness_nonlinear_plane_strain!(K_nonlinear, dh, cellvalues, input::InputStruct, ue)
    n_basefuncs = getnbasefunctions(cellvalues)
    ke_n = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K_nonlinear)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(ke_n, 0.0)
        local_stiffness_nonlinear_plane_strain!(ke_n, cell, cellvalues, input, ue)
        assemble!(assembler, celldofs(cell), ke_n)
    end
    return K_nonlinear
end
########################### end of global_siffness_nonlinear! ##################################
"""
assemble_internal_force_plane_stain!(f_int::Vector, cell::CellCache, cv::CellValues, input::InputStruct, ue_global::Vector)

Internal force for nonlinear part
"""
function assemble_internal_force_plane_stain!(f_int::Vector, cell::CellCache, cv::CellValues, input::InputStruct, ue_global::Vector)
    reinit!(cv, cell)
    fill!(f_int, 0.0)
    ndofs = getnbasefunctions(cv)
    material = input.material
    cell_dofs = celldofs(cell)
    ue_local = view(ue_global, cell_dofs)
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇u2d = function_gradient(cv, qp, ue_local)
        F = Tensor{2,3,Float64}((
            1.0 + ∇u2d[1, 1], ∇u2d[1, 2], 0.0,
            ∇u2d[2, 1], 1.0 + ∇u2d[2, 2], 0.0,
            0.0, 0.0, 1.0
        ))
        C = tdot(F)
        S, _ = material(C)
        P = F ⋅ S
        for i in 1:ndofs
            ∇δui2d = shape_gradient(cv, qp, i)
            ∇δui = Tensor{2,3,Float64}((
                ∇δui2d[1, 1], ∇δui2d[1, 2], 0.0,
                ∇δui2d[2, 1], ∇δui2d[2, 2], 0.0,
                0.0, 0.0, 0.0
            ))
            f_int[i] += (∇δui ⊡ P) * dΩ
        end
    end
    return
end
######################################### end of assemble_internal_force! ######################################
"""
assemble_global_internal_force_plane_strain!(f_int_global, dh, cellvalues, input::InputStruct, ue)

Assemble internal force stiffness
"""
function assemble_global_internal_force_plane_strain!(f_int_global, dh, cellvalues, input::InputStruct, ue)
    n_basefuncs = getnbasefunctions(cellvalues)
    f_int_local = zeros(n_basefuncs)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(f_int_local, 0.0)
        assemble_internal_force_plane_stain!(f_int_local, cell, cellvalues, input, ue)
        assemble!(f_int_global, celldofs(cell), f_int_local)
    end
    return f_int_global
end
######################### end of assemble_global_internal_force! ###########################
"""
assemble_traction_forces!(f_ext, dh, facetset, facetvalues, prescribed_traction)

Assemble for traction force (surface integral)
"""
# function assemble_traction_forces!(f_ext, dh, facetset, facetvalues, prescribed_traction)
#     fe_ext = zeros(getnbasefunctions(facetvalues))
#     for facet in FacetIterator(dh, facetset)
#         reinit!(facetvalues, facet)
#         fill!(fe_ext, 0.0)
#         for qp in 1:getnquadpoints(facetvalues)
#             # Calculate the global coordinate of the quadrature point.
#             tₚ = prescribed_traction
#             dΓ = getdetJdV(facetvalues, qp)
#             for i in 1:getnbasefunctions(facetvalues)
#                 Nᵢ = shape_value(facetvalues, qp, i)
#                 fe_ext[i] += tₚ ⋅ Nᵢ * dΓ
#             end
#         end
#         assemble!(f_ext, celldofs(facet), fe_ext)
#     end
#     return f_ext
# end

function assemble_traction_forces!(
    f_ext,
    dh,
    facetsets::Vector,
    facetvalues,
    tractions::Dict{Int, Vector}
)
    fe_ext = zeros(getnbasefunctions(facetvalues))
    for (idx, facetset) in enumerate(facetsets)
        t_traction = tractions[idx]
        for facet in FacetIterator(dh, facetset)
            reinit!(facetvalues, facet)
            fill!(fe_ext, 0.0)
            for qp in 1:getnquadpoints(facetvalues)
                dΓ = getdetJdV(facetvalues, qp)
                for i in 1:getnbasefunctions(facetvalues)
                    Nᵢ = shape_value(facetvalues, qp, i)
                    fe_ext[i] += t_traction ⋅ Nᵢ * dΓ
                end
            end
            assemble!(f_ext, celldofs(facet), fe_ext)
        end
    end
    return f_ext
end

###### for multiple load cases

################# end of assemble_traction_forces! ################################
"""
solve_lambda3(F2d, input::InputStruct; tol=1e-10, maxit=25)
--- Plane stress: Solve for λ3 given in-plane F2d so that S33 = 0 ---
"""

function solve_lambda3(F2d, input::InputStruct; tol=1e-10, maxit=25)
    material = input.material
    λ3 = 1.0
    for k in 1:maxit
        F = Tensor{2,3,Float64}((
            F2d[1,1], F2d[1,2], 0.0,
            F2d[2,1], F2d[2,2], 0.0,
            0.0,      0.0,      λ3
        ))
        C = tdot(F)
        S, _ = material(C)
        S33 = S[3,3]
        δ = 1e-8
        λ3p = λ3 + δ
        Fp = Tensor{2,3,Float64}((
            F2d[1,1], F2d[1,2], 0.0,
            F2d[2,1], F2d[2,2], 0.0,
            0.0,      0.0,      λ3p
        ))
        Cp = tdot(Fp)
        Sp, _ = material(Cp)
        S33p = Sp[3,3]
        dS33dλ3 = (S33p - S33)/δ
        λ3new = λ3 - S33/dS33dλ3
        if abs(λ3new - λ3) < tol
            return λ3new
        end
        λ3 = λ3new
    end
    error("λ3 for plane stress did not converge")
end
############################# end of solve_lambda3 ######################
"""
local_stiffness_nonlinear_plane_stress!(ke_n, cell,cv, input::InputStruct, ue_global)

Nonlinear assembly (plane stress support) ---
"""
function local_stiffness_nonlinear_plane_stress!(ke_n, cell,cv, input::InputStruct, ue_global)
    reinit!(cv, cell)
    fill!(ke_n, 0.0)
    ndofs = getnbasefunctions(cv)
    cell_dofs = celldofs(cell)
    ue_local = view(ue_global, cell_dofs)
    material = input.material
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇u2d = function_gradient(cv, qp, ue_local)
        # Compose F2d (2x2)
        F2d = [1.0 + ∇u2d[1, 1]  ∇u2d[1, 2];
               ∇u2d[2, 1]        1.0 + ∇u2d[2, 2]]
        # Plane stress: solve for λ3
        λ3 = solve_lambda3(F2d, input)
        # Build full F
        F = Tensor{2,3,Float64}((
            F2d[1,1], F2d[1,2], 0.0,
            F2d[2,1], F2d[2,2], 0.0,
            0.0,      0.0,      λ3
        ))
        C = tdot(F)
        S, ∂S∂C = material(C)
        I = one(S)
        ∂P∂F = otimesu(I, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(F', I)
        for i in 1:ndofs
            ∇δui2d = shape_gradient(cv, qp, i)
            ∇δui = Tensor{2,3,Float64}((
                ∇δui2d[1,1], ∇δui2d[1,2], 0.0,
                ∇δui2d[2,1], ∇δui2d[2,2], 0.0,
                0.0,         0.0,         0.0
            ))
            ∇δui_∂P∂F = ∇δui ⊡ ∂P∂F
            for j in 1:ndofs
                ∇δuj2d = shape_gradient(cv, qp, j)
                ∇δuj = Tensor{2,3,Float64}((
                    ∇δuj2d[1,1], ∇δuj2d[1,2], 0.0,
                    ∇δuj2d[2,1], ∇δuj2d[2,2], 0.0,
                    0.0,         0.0,         0.0
                ))
                ke_n[i, j] += (∇δui_∂P∂F ⊡ ∇δuj) * dΩ
            end
        end
    end
    return
end
########################### end of local_stiffness_nonlinear_plane_stress! ##################
"""
global_siffness_nonlinear_plane_stress!(K_nonlinear, dh, cellvalues, input::InputStruct, ue)
Assemble local stiffness nonlinear part
"""
function global_siffness_nonlinear_plane_stress!(K_nonlinear, dh, cellvalues, input::InputStruct, ue)
    n_basefuncs = getnbasefunctions(cellvalues)
    ke_n = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K_nonlinear)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(ke_n, 0.0)
        local_stiffness_nonlinear_plane_stress!(ke_n, cell, cellvalues, input, ue)
        assemble!(assembler, celldofs(cell), ke_n)
    end
    return K_nonlinear
end
######################### end of global_siffness_nonlinear_plane_stress! ################
""""
assemble_internal_force_plane_stress!(f_int::Vector, cell::CellCache, cv::CellValues, input::InputStruct, ue_global::Vector)
# --- Internal force assembly (plane stress) ---
"""
function assemble_internal_force_plane_stress!(f_int::Vector, cell::CellCache, cv::CellValues, input::InputStruct, ue_global::Vector)
    reinit!(cv, cell)
    fill!(f_int, 0.0)
    ndofs = getnbasefunctions(cv)
    cell_dofs = celldofs(cell)
    ue_local = view(ue_global, cell_dofs)
    material = input.material
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇u2d = function_gradient(cv, qp, ue_local)
        F2d = [1.0 + ∇u2d[1, 1]  ∇u2d[1, 2];
               ∇u2d[2, 1]        1.0 + ∇u2d[2, 2]]
        λ3 = solve_lambda3(F2d, input)
        F = Tensor{2,3,Float64}((
            F2d[1,1], F2d[1,2], 0.0,
            F2d[2,1], F2d[2,2], 0.0,
            0.0,      0.0,      λ3
        ))
        C = tdot(F)
        S, ∂S∂C = material(C)
        P = F ⋅ S
        for i in 1:ndofs
            ∇δui2d = shape_gradient(cv, qp, i)
            ∇δui = Tensor{2,3,Float64}((
                ∇δui2d[1,1], ∇δui2d[1,2], 0.0,
                ∇δui2d[2,1], ∇δui2d[2,2], 0.0,
                0.0,         0.0,         0.0
            ))
            f_int[i] += (∇δui ⊡ P) * dΩ
        end
    end
    return
end
############################ End of assemble_internal_force_plane_stress! #############
"""
assemble_global_internal_force_plane_stress!(f_int_global, dh, cellvalues, input::InputStruct, ue)
Assemble global nonlinear part plane stress
"""
function assemble_global_internal_force_plane_stress!(f_int_global, dh, cellvalues, input::InputStruct, ue)
    n_basefuncs = getnbasefunctions(cellvalues)
    f_int_local = zeros(n_basefuncs)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(f_int_local, 0.0)
        assemble_internal_force_plane_stress!(f_int_local, cell, cellvalues, input, ue)
        assemble!(f_int_global, celldofs(cell), f_int_local)
    end
    return f_int_global
end
############################### end of assemble_global_internal_force_plane_stress! ################
####### 3D cases, nonlinear part
"""
local_stiffness_nonlinear_3D!(ke_n, cell, cv, input::InputStruct, ue_global)
Local stiffness 3D cases
"""
function local_stiffness_nonlinear_3D!(ke_n, cell, cv, input::InputStruct, ue_global)
    reinit!(cv, cell)
    fill!(ke_n, 0.0)

    ndofs = getnbasefunctions(cv)
    cell_dofs = celldofs(cell)
    ue_local = view(ue_global, cell_dofs)
    material = input.material

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)

        # Compute deformation gradient F and right Cauchy-Green tensor C
        ∇u = function_gradient(cv, qp, ue_local)  # Tensor{2,3}
        F = one(∇u) + ∇u
        C = tdot(F)

        # Get stress and tangent from material model
        S, ∂S∂C = material(C)

        I = one(S)
        ∂P∂F = otimesu(I, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(transpose(F), I)

        for i in 1:ndofs
            ∇δui = shape_gradient(cv, qp, i)  # Tensor{2,3}
            ∇δui_∂P∂F = ∇δui ⊡ ∂P∂F

            for j in 1:ndofs
                ∇δuj = shape_gradient(cv, qp, j)  # Tensor{2,3}
                ke_n[i, j] += (∇δui_∂P∂F ⊡ ∇δuj) * dΩ
            end
        end
    end

    return
end
########################### end of local_stiffness_nonlinear_3D! ####################
"""
global_siffness_nonlinear_3D!(K_nonlinear, dh, cellvalues, input::InputStruct, ue)

"""
function global_siffness_nonlinear_3D!(K_nonlinear, dh, cellvalues, input::InputStruct, ue)
    n_basefuncs = getnbasefunctions(cellvalues)
    ke_n = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K_nonlinear)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(ke_n, 0.0)
        local_stiffness_nonlinear_3D!(ke_n, cell, cellvalues, input, ue)
        assemble!(assembler, celldofs(cell), ke_n)
    end
    return K_nonlinear
end
################################## end of global_siffness_nonlinear_3D! ################
"""
assemble_internal_force_3D!(f_int::Vector, cell::CellCache, cv::CellValues, input::InputStruct, ue_global::Vector)

"""
function assemble_internal_force_3D!(f_int::Vector, cell::CellCache, cv::CellValues, input::InputStruct, ue_global::Vector)
    reinit!(cv, cell)
    fill!(f_int, 0.0)

    ndofs = getnbasefunctions(cv)
    material = input.material
    cell_dofs = celldofs(cell)
    ue_local = view(ue_global, cell_dofs)

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)

        # Deformation gradient F
        ∇u = function_gradient(cv, qp, ue_local)  # 3×3 tensor
        F = one(∇u) + ∇u
        C = tdot(F)

        # Material response (stress)
        S, _ = material(C)
        P = F ⋅ S   # 1st Piola-Kirchhoff stress

        # Internal force vector
        for i in 1:ndofs
            ∇δui = shape_gradient(cv, qp, i)  # 3×3 tensor
            f_int[i] += (∇δui ⊡ P) * dΩ
        end
    end

    return
end
################################# end of assemble_internal_force_3D! ###################
"""
assemble_global_internal_force_3D!(f_int_global, dh, cellvalues, input::InputStruct, ue)
"""
function assemble_global_internal_force_3D!(f_int_global, dh, cellvalues, input::InputStruct, ue)
    n_basefuncs = getnbasefunctions(cellvalues)
    f_int_local = zeros(n_basefuncs)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(f_int_local, 0.0)
        assemble_internal_force_3D!(f_int_local, cell, cellvalues, input, ue)
        assemble!(f_int_global, celldofs(cell), f_int_local)
    end
    return f_int_global
end
########################### end of assemble_global_internal_force_3D!
#######################################################

struct FESolver
    u::Vector{Float64}
end

"""
run_plane_strain(input::InputStruct)
"""
function run_plane_strain(input::InputStruct)
    #grid = input.grid
    cv = input.cell_values
    fv = input.facet_values
    dh = input.dh
    ch = input.ch
    ΓN = input.facetsets
    tol = input.tol
    n_load_steps = input.n_load_steps
    traction = input.tractions
    filename = input.filename
    output_dir = input.output_dir
    n_iter_NR = input.n_iter_NR
    # --- Initialization ---
    E = input.E
    ν = input.ν
    ElaTensor = elasticity_tensor_plane_strain(E, ν)

    U = zeros(ndofs(dh))           # Final total displacement
    Uinit = zeros(ndofs(dh))       # Running displacement for current load step
    Residual = zeros(ndofs(dh))
    conv_flag = false

    for O in 1:n_load_steps
        println("\n🔷 Load step $O")

        # traction_load = traction /n_load_steps
        traction_load = Dict{Int, Vector}(k => v ./ n_load_steps for (k, v) in traction)
        f_ext = zeros(ndofs(dh))
        assemble_traction_forces!(f_ext, dh, ΓN, fv, traction_load)

        K_linear = allocate_matrix(dh)
        global_stiffness_linear!(K_linear, dh, cv, ElaTensor)
        update!(ch,O )
        apply!(K_linear, f_ext, ch)
        # Initial linear guess
        dU = K_linear \ f_ext
        Uinit .= Uinit .+ dU

        # --- Newton-Raphson Iteration ---
        for L in 1:n_iter_NR
            # Assemble nonlinear stiffness matrix
            K_nonlinear = allocate_matrix(dh)
            global_siffness_nonlinear_plane_strain!(K_nonlinear, dh, cv, input, Uinit)

            f_int_global = zeros(ndofs(dh))
            assemble_global_internal_force_plane_strain!(f_int_global, dh, cv, input, Uinit)

            Residual .= f_ext - f_int_global
            update!(ch,O)
            # Apply BCs to residual and stiffness matrix
            apply!(K_nonlinear, Residual, ch)
            # Solve for update
            deltaU = zeros(ndofs(dh))
            deltaU .= K_nonlinear \ Residual

            # Update displacement
            Uinit .= Uinit .+ deltaU

            # Check convergence
            conv = norm(Residual) / (1 + norm(f_ext))
            #println("  🔁 Newton iter $L: residual = $(round(conv, sigdigits=4))")
            if conv < tol
                conv_flag = true
                break
            else
                conv_flag = false
            end
        end

        # Check if convergence was achieved
        if !conv_flag
            println("❌ CONVERGENCE DID NOT REACH FOR STEP: $O")
        else
            println("✅ CONVERGENCE IS REACHED FOR STEP: $O")
            U .= U .+ Uinit
        end
        VTKGridFile(joinpath(output_dir,filename), dh) do vtk
            write_solution(vtk, dh, U)
        end
    end

    return FESolver(U)
end
######################## end of run_plane_strain
function run_plane_stress(input::InputStruct)
    #grid = input.grid
    cv = input.cell_values
    fv = input.facet_values
    dh = input.dh
    ch = input.ch
    ΓN = input.facetsets
    tol = input.tol
    n_load_steps = input.n_load_steps
    traction = input.tractions
    filename = input.filename
    output_dir = input.output_dir
    n_iter_NR = input.n_iter_NR
    # --- Initialization ---
    E = input.E
    ν = input.ν
    ElaTensor = elasticity_tensor_plane_stress(E, ν)

    U = zeros(ndofs(dh))           # Final total displacement
    Uinit = zeros(ndofs(dh))       # Running displacement for current load step
    Residual = zeros(ndofs(dh))
    conv_flag = false

    for O in 1:n_load_steps
        println("\n🔷 Load step $O")

        # traction_load = traction /n_load_steps
        traction_load = Dict{Int, Vector}(k => v ./ n_load_steps for (k, v) in traction)
        f_ext = zeros(ndofs(dh))
        assemble_traction_forces!(f_ext, dh, ΓN, fv, traction_load)

        K_linear = allocate_matrix(dh)
        global_stiffness_linear!(K_linear, dh, cv, ElaTensor)
        update!(ch,O )
        apply!(K_linear, f_ext, ch)
        # Initial linear guess
        dU = K_linear \ f_ext
        Uinit .= Uinit .+ dU

        # --- Newton-Raphson Iteration ---
        for L in 1:n_iter_NR
            # Assemble nonlinear stiffness matrix
            K_nonlinear = allocate_matrix(dh)
            global_siffness_nonlinear_plane_stress!(K_nonlinear, dh, cv, input, Uinit)

            f_int_global = zeros(ndofs(dh))
            assemble_global_internal_force_plane_stress!(f_int_global, dh, cv, input, Uinit)

            Residual .= f_ext - f_int_global
            update!(ch,O)
            # Apply BCs to residual and stiffness matrix
            apply!(K_nonlinear, Residual, ch)
            # Solve for update
            deltaU = zeros(ndofs(dh))
            deltaU .= K_nonlinear \ Residual

            # Update displacement
            Uinit .= Uinit .+ deltaU

            # Check convergence
            conv = norm(Residual) / (1 + norm(f_ext))
            #println("  🔁 Newton iter $L: residual = $(round(conv, sigdigits=4))")
            if conv < tol
                conv_flag = true
                break
            else
                conv_flag = false
            end
        end

        # Check if convergence was achieved
        if !conv_flag
            println("❌ CONVERGENCE DID NOT REACH FOR STEP: $O")
        else
            println("✅ CONVERGENCE IS REACHED FOR STEP: $O")
            U .= U .+ Uinit
        end
        VTKGridFile(joinpath(output_dir,filename), dh) do vtk
            write_solution(vtk, dh, U)
        end
    end

    return FESolver(U)
end
################ end of run_plane_stress
function run_3d(input::InputStruct)
    #grid = input.grid
    cv = input.cell_values
    fv = input.facet_values
    dh = input.dh
    ch = input.ch
    ΓN = input.facetsets
    tol = input.tol
    n_load_steps = input.n_load_steps
    traction = input.tractions
    filename = input.filename
    output_dir = input.output_dir
    n_iter_NR = input.n_iter_NR
    # --- Initialization ---
    E = input.E
    ν = input.ν
    ElaTensor = elasticity_tensor_3D(E, ν)

    U = zeros(ndofs(dh))           # Final total displacement
    Uinit = zeros(ndofs(dh))       # Running displacement for current load step
    Residual = zeros(ndofs(dh))
    conv_flag = false

    for O in 1:n_load_steps
        println("\n🔷 Load step $O")

        # traction_load = traction /n_load_steps
        traction_load = Dict{Int, Vector}(k => v ./ n_load_steps for (k, v) in traction)
        f_ext = zeros(ndofs(dh))
        assemble_traction_forces!(f_ext, dh, ΓN, fv, traction_load)

        K_linear = allocate_matrix(dh)
        global_stiffness_linear!(K_linear, dh, cv, ElaTensor)
        update!(ch,O )
        apply!(K_linear, f_ext, ch)
        # Initial linear guess
        dU = K_linear \ f_ext
        Uinit .= Uinit .+ dU

        # --- Newton-Raphson Iteration ---
        for L in 1:n_iter_NR
            # Assemble nonlinear stiffness matrix
            K_nonlinear = allocate_matrix(dh)
            global_siffness_nonlinear_3D!(K_nonlinear, dh, cv, input, Uinit)

            f_int_global = zeros(ndofs(dh))
            assemble_global_internal_force_3D!(f_int_global, dh, cv, input, Uinit)

            Residual .= f_ext - f_int_global
            update!(ch,O)
            # Apply BCs to residual and stiffness matrix
            apply!(K_nonlinear, Residual, ch)
            # Solve for update
            deltaU = zeros(ndofs(dh))
            deltaU .= K_nonlinear \ Residual

            # Update displacement
            Uinit .= Uinit .+ deltaU

            # Check convergence
            conv = norm(Residual) / (1 + norm(f_ext))
            #println("  🔁 Newton iter $L: residual = $(round(conv, sigdigits=4))")
            if conv < tol
                conv_flag = true
                break
            else
                conv_flag = false
            end
        end

        # Check if convergence was achieved
        if !conv_flag
            println("❌ CONVERGENCE DID NOT REACH FOR STEP: $O")
        else
            println("✅ CONVERGENCE IS REACHED FOR STEP: $O")
            U .= U .+ Uinit
        end
        VTKGridFile(joinpath(output_dir,filename), dh) do vtk
            write_solution(vtk, dh, U)
        end
    end

    return FESolver(U)
end

############## end of run_3d
"""
run_fem(input::InputStruct)
"""
function run_fem(input::InputStruct)
    if input.model_type == :plane_stress
        run_plane_stress(input)
    elseif input.model_type == :plane_strain
        run_plane_strain(input)
    elseif input.model_type == :threeD
        run_3d(input)
    else
        error("Unknown model_type: $(input.model_type)")
    end
end
