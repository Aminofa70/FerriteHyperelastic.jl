using LinearAlgebra, Printf
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
local_stiffness_nonlinear!(ke_n, cell, cv, input::InputStruct, ue_global)
Local Stiffness Matrix for nonlinear part: plane strain
"""
function local_stiffness_nonlinear!(ke_n, cell, cv, input::InputStruct, ue_global)
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
global_siffness_nonlinear!(K_nonlinear, dh, cellvalues, input::InputStruct, ue)

"""
function global_siffness_nonlinear!(K_nonlinear, dh, cellvalues, input::InputStruct, ue)
    n_basefuncs = getnbasefunctions(cellvalues)
    ke_n = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K_nonlinear)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(ke_n, 0.0)
        local_stiffness_nonlinear!(ke_n, cell, cellvalues, input, ue)
        assemble!(assembler, celldofs(cell), ke_n)
    end
    return K_nonlinear
end
########################### end of global_siffness_nonlinear! ##################################
"""
assemble_internal_force!(f_int::Vector, cell::CellCache, cv::CellValues, input::InputStruct, ue_global::Vector)

Internal force for nonlinear part
"""
function assemble_internal_force!(f_int::Vector, cell::CellCache, cv::CellValues, input::InputStruct, ue_global::Vector)
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
assemble_global_internal_force!(f_int_global, dh, cellvalues, input::InputStruct, ue)

Assemble internal force stiffness
"""
function assemble_global_internal_force!(f_int_global, dh, cellvalues, input::InputStruct, ue)
    n_basefuncs = getnbasefunctions(cellvalues)
    f_int_local = zeros(n_basefuncs)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(f_int_local, 0.0)
        assemble_internal_force!(f_int_local, cell, cellvalues, input, ue)
        assemble!(f_int_global, celldofs(cell), f_int_local)
    end
    return f_int_global
end
######################### end of assemble_global_internal_force! ###########################
"""
assemble_traction_forces!(f_ext, dh, facetset, facetvalues, prescribed_traction)

Assemble for traction force (surface integral)
"""
function FerriteHyperelastic.assemble_traction_forces!(f_ext, dh, facetset, facetvalues, prescribed_traction)
    fe_ext = zeros(getnbasefunctions(facetvalues))
    for facet in FacetIterator(dh, facetset)
        reinit!(facetvalues, facet)
        fill!(fe_ext, 0.0)
        for qp in 1:getnquadpoints(facetvalues)
            # Calculate the global coordinate of the quadrature point.
            tₚ = prescribed_traction
            dΓ = getdetJdV(facetvalues, qp)
            for i in 1:getnbasefunctions(facetvalues)
                Nᵢ = shape_value(facetvalues, qp, i)
                fe_ext[i] += tₚ ⋅ Nᵢ * dΓ
            end
        end
        assemble!(f_ext, celldofs(facet), fe_ext)
    end
    return f_ext
end
################# end of assemble_traction_forces! ################################
struct FESolver
    u::Vector{Float64}
end

function run_fem(input::InputStruct)
    grid = input.grid
    cv = input.cell_values
    fv = input.facet_values
    dh = input.dh
    ch = input.ch
    ΓN = par.ΓN
    tol = input.tol
    n_load_steps = input.n_load_steps
    traction = input.traction
    filename = input.filename
    output_dir = input.output_dir
    n_iter_NR = input.n_iter_NR
    # --- Initialization ---


    U = zeros(ndofs(dh))           # Final total displacement
    Uinit = zeros(ndofs(dh))       # Running displacement for current load step
    Residual = zeros(ndofs(dh))
    conv_flag = false

    for O in 1:n_load_steps
        println("\n🔷 Load step $O")

        traction_load = traction /n_load_steps
        f_ext = zeros(ndofs(dh))
        assemble_external_forces!(f_ext, dh, ΓN, fv, traction_load)

        K_linear = allocate_matrix(dh)
        global_stiffness_linear!(K_linear, dh, cv, input)
        update!(ch,O )
        apply!(K_linear, f_ext, ch)
        # Initial linear guess
        dU = K_linear \ f_ext
        Uinit .= Uinit .+ dU

        # --- Newton-Raphson Iteration ---
        for L in 1:n_iter_NR
            # Assemble nonlinear stiffness matrix
            K_nonlinear = allocate_matrix(dh)
            global_siffness_nonlinear!(K_nonlinear, dh, cellvalues, input, ue)

            f_int_global = zeros(ndofs(dh))
            assemble_global_internal_force!(f_int_global, dh, cv, input, Uinit)

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