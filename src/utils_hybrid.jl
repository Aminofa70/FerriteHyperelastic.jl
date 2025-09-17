function assemble_element_plane_strain!(Ke_n, fe_int, cell, cellvalues_u, cellvalues_p, input::InputStruct, ue, pe)

    # Reinitialize cell values, and reset output arrays
    ublock, pblock = 1, 2
    Ferrite.reinit!(cellvalues_u, cell)
    Ferrite.reinit!(cellvalues_p, cell)
    fill!(Ke_n, 0.0)
    fill!(fe_int, 0.0)

    n_basefuncs_u = getnbasefunctions(cellvalues_u)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)

    material = input.material

    for qp in 1:getnquadpoints(cellvalues_u)
        dΩ = getdetJdV(cellvalues_u, qp)
        # Compute deformation gradient F

        p = function_value(cellvalues_p, qp, pe)
        ∇u2d = function_gradient(cellvalues_u, qp, ue)
        F = Tensor{2,3,Float64}((
            1.0 + ∇u2d[1, 1], ∇u2d[1, 2], 0.0,
            ∇u2d[2, 1], 1.0 + ∇u2d[2, 2], 0.0,
            0.0, 0.0, 1.0
        ))

        # Compute first Piola-Kirchhoff stress and tangent modulus
        ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂p, ∂²Ψ∂p², ∂²Ψ∂F∂p = material(F, p)

        # Loop over the `u`-test functions to calculate the `u`-`u` and `u`-`p` blocks
        for i in 1:n_basefuncs_u
            # gradient of the test function

            ∇δui2d = shape_gradient(cellvalues_u, qp, i)
            ∇δui = Tensor{2,3,Float64}((
                ∇δui2d[1, 1], ∇δui2d[1, 2], 0.0,
                ∇δui2d[2, 1], ∇δui2d[2, 2], 0.0,
                0.0, 0.0, 0.0
            ))

            # Add contribution to the residual from this test function
            fe_int[BlockIndex((ublock), (i))] += (∇δui ⊡ ∂Ψ∂F) * dΩ
            ∇δui∂S∂F = ∇δui ⊡ ∂²Ψ∂F²



            for j in 1:n_basefuncs_u

                ∇δuj2d = shape_gradient(cellvalues_u, qp, j)
                ∇δuj = Tensor{2,3,Float64}((
                    ∇δuj2d[1, 1], ∇δuj2d[1, 2], 0.0,
                    ∇δuj2d[2, 1], ∇δuj2d[2, 2], 0.0,
                    0.0, 0.0, 0.0
                ))
                # Add contribution to the tangent
                Ke_n[BlockIndex((ublock, ublock), (i, j))] += (∇δui∂S∂F ⊡ ∇δuj) * dΩ

            end
            # Loop over the `p`-test functions
            for j in 1:n_basefuncs_p
                δp = shape_value(cellvalues_p, qp, j)
                # Add contribution to the tangent
                Ke_n[BlockIndex((ublock, pblock), (i, j))] += (∂²Ψ∂F∂p ⊡ ∇δui) * δp * dΩ
            end
        end
        # Loop over the `p`-test functions to calculate the `p-`u` and `p`-`p` blocks
        for i in 1:n_basefuncs_p
            δp = shape_value(cellvalues_p, qp, i)
            fe_int[BlockIndex((pblock), (i))] += (δp * ∂Ψ∂p) * dΩ

            for j in 1:n_basefuncs_u



                ∇δuj2d = shape_gradient(cellvalues_u, qp, j)
                ∇δuj = Tensor{2,3,Float64}((
                    ∇δuj2d[1, 1], ∇δuj2d[1, 2], 0.0,
                    ∇δuj2d[2, 1], ∇δuj2d[2, 2], 0.0,
                    0.0, 0.0, 0.0
                ))
                Ke_n[BlockIndex((pblock, ublock), (i, j))] += ∇δuj ⊡ ∂²Ψ∂F∂p * δp * dΩ


            end
            for j in 1:n_basefuncs_p
                δp = shape_value(cellvalues_p, qp, j)
                Ke_n[BlockIndex((pblock, pblock), (i, j))] += δp * ∂²Ψ∂p² * δp * dΩ
            end
        end
    end
    return
end;
####################################################################
function assemble_element_threeD!(Ke_n, fe_int, cell, cellvalues_u, cellvalues_p, input::InputStruct, ue, pe)

    # Reinitialize cell values, and reset output arrays
    ublock, pblock = 1, 2
    Ferrite.reinit!(cellvalues_u, cell)
    Ferrite.reinit!(cellvalues_p, cell)
    fill!(Ke_n, 0.0)
    fill!(fe_int, 0.0)

    n_basefuncs_u = getnbasefunctions(cellvalues_u)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)

    material = input.material

    for qp in 1:getnquadpoints(cellvalues_u)
        dΩ = getdetJdV(cellvalues_u, qp)
        # Compute deformation gradient F

        p = function_value(cellvalues_p, qp, pe)
        ∇u = function_gradient(cellvalues_u, qp, ue)
        F = one(∇u) + ∇u

        # Compute first Piola-Kirchhoff stress and tangent modulus
        ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂p, ∂²Ψ∂p², ∂²Ψ∂F∂p = material(F, p)

        # Loop over the `u`-test functions to calculate the `u`-`u` and `u`-`p` blocks
        for i in 1:n_basefuncs_u
            # gradient of the test function

            ∇δui = shape_gradient(cellvalues_u, qp, i)
        
            # Add contribution to the residual from this test function
            fe_int[BlockIndex((ublock), (i))] += (∇δui ⊡ ∂Ψ∂F) * dΩ
            ∇δui∂S∂F = ∇δui ⊡ ∂²Ψ∂F²

            for j in 1:n_basefuncs_u
                ∇δuj = shape_gradient(cellvalues_u, qp, j)
                # Add contribution to the tangent
                Ke_n[BlockIndex((ublock, ublock), (i, j))] += (∇δui∂S∂F ⊡ ∇δuj) * dΩ
            end
            # Loop over the `p`-test functions
            for j in 1:n_basefuncs_p
                δp = shape_value(cellvalues_p, qp, j)
                # Add contribution to the tangent
                Ke_n[BlockIndex((ublock, pblock), (i, j))] += (∂²Ψ∂F∂p ⊡ ∇δui) * δp * dΩ
            end
        end
        # Loop over the `p`-test functions to calculate the `p-`u` and `p`-`p` blocks
        for i in 1:n_basefuncs_p
            δp = shape_value(cellvalues_p, qp, i)
            fe_int[BlockIndex((pblock), (i))] += (δp * ∂Ψ∂p) * dΩ

            for j in 1:n_basefuncs_u
                ∇δuj = shape_gradient(cellvalues_u, qp, j)
                Ke_n[BlockIndex((pblock, ublock), (i, j))] += ∇δuj ⊡ ∂²Ψ∂F∂p * δp * dΩ
            end
            for j in 1:n_basefuncs_p
                δp = shape_value(cellvalues_p, qp, j)
                Ke_n[BlockIndex((pblock, pblock), (i, j))] += δp * ∂²Ψ∂p² * δp * dΩ
            end
        end
    end
    return
end;
####################################################################
function assemble_global!(K_n::SparseMatrixCSC, f_int, cellvalues_u::CellValues,
    cellvalues_p::CellValues, dh::DofHandler, input::InputStruct, w
)
    nu = getnbasefunctions(cellvalues_u)
    np = getnbasefunctions(cellvalues_p)

    # start_assemble resets K and f
    fe_int = BlockedArray(zeros(nu + np), [nu, np]) # local force vector
    ke_n = BlockedArray(zeros(nu + np, nu + np), [nu, np], [nu, np]) # local stiffness matrix

    assembler = start_assemble(K_n, f_int)
    # Loop over all cells in the grid
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        global_dofsu = global_dofs[1:nu] # first nu dofs are displacement
        global_dofsp = global_dofs[(nu+1):end] # last np dofs are pressure
        @assert size(global_dofs, 1) == nu + np # sanity check
        ue = w[global_dofsu] # displacement dofs for the current cell
        pe = w[global_dofsp] # pressure dofs for the current cell
        if input.model_type == :plane_strain
            assemble_element_plane_strain!(ke_n, fe_int, cell, cellvalues_u, cellvalues_p, input, ue, pe)
        elseif input.model_type == :threeD
            assemble_element_threeD!(ke_n, fe_int, cell, cellvalues_u, cellvalues_p, input, ue, pe)
        else
            error("the model type is not correct")
        end
        assemble!(assembler, global_dofs, ke_n, fe_int)
    end
    return
end;
####################################################################
function assemble_traction_forces_twoD_hybrid!(F_ext, dh, facetsets, facetvalues, tractions, u, cellvalues_u, cellvalues_p)
    nu = getnbasefunctions(cellvalues_u) # Number of displacement DOFs
    np = getnbasefunctions(cellvalues_p) # Number of pressure DOFs
    fe_ext = BlockedArray(zeros(nu + np), [nu, np]) # Local force vector with u and p blocks
    for (idx, facetset) in enumerate(facetsets)
        t_traction = Ferrite.Vec{2,Float64}(tractions[idx]) # Ensure traction is a 2D vector
        for facetidx in facetset
            cellid = facetidx[1]
            facet_idx = facetidx[2]
            cell = getcells(dh.grid)[cellid]
            cell_coords = getcoordinates(dh.grid, cellid)
            updated_coords = copy(cell_coords)
            cell_dofs = celldofs(dh, cellid)
            u_dof = dof_range(dh, :u) # Displacement DOF indices
            for i in eachindex(cell_coords)
                current_coords = cell_coords[i]
                dof_x = cell_dofs[u_dof[2*i-1]] # x-displacement DOF
                dof_y = cell_dofs[u_dof[2*i]]   # y-displacement DOF
                new_coords = (current_coords[1] + u[dof_x], current_coords[2] + u[dof_y])
                updated_coords[i] = Ferrite.Vec{2,Float64}(new_coords)
            end
            Ferrite.reinit!(facetvalues, cell, updated_coords, facet_idx)
            fill!(fe_ext, 0.0) # Reset local force vector
            for qp in 1:getnquadpoints(facetvalues)
                dΓ = getdetJdV(facetvalues, qp)
                for i in 1:getnbasefunctions(facetvalues)
                    Nᵢ = shape_value(facetvalues, qp, i) # Shape function (vector for displacement)
                    # Add traction contribution to displacement block
                    fe_ext[BlockIndex((1), (i))] += (t_traction ⋅ Nᵢ) * dΓ
                end
            end
            # Assemble into global force vector using all DOFs (u and p)
            assemble!(F_ext, cell_dofs, fe_ext)
        end
    end
    return F_ext
end
####################################################################
function assemble_traction_forces_threeD_hybrid!(F_ext, dh, facetsets, facetvalues, tractions, u, cellvalues_u, cellvalues_p)
    nu = getnbasefunctions(cellvalues_u) # Number of displacement DOFs
    np = getnbasefunctions(cellvalues_p) # Number of pressure DOFs
    fe_ext = BlockedArray(zeros(nu + np), [nu, np]) # Local force vector with u and p blocks
    for (idx, facetset) in enumerate(facetsets)
        t_traction = Ferrite.Vec{3,Float64}(tractions[idx]) # Ensure traction is a 2D vector
        for facetidx in facetset
            cellid = facetidx[1]
            facet_idx = facetidx[2]
            cell = getcells(dh.grid)[cellid]
            cell_coords = getcoordinates(dh.grid, cellid)
            updated_coords = copy(cell_coords)
            cell_dofs = celldofs(dh, cellid)
            u_dof = dof_range(dh, :u) # Displacement DOF indices
            for i in eachindex(cell_coords)
                current_coords = cell_coords[i]
                # Get DOFs for the i-th node (x, y, z displacements)
                dof_x = cell_dofs[3*i-2]  # x-displacement DOF
                dof_y = cell_dofs[3*i-1]  # y-displacement DOF
                dof_z = cell_dofs[3*i]    # z-displacement DOF
                # Create new coordinates
                new_coords = (
                    current_coords[1] + u[dof_x],
                    current_coords[2] + u[dof_y],
                    current_coords[3] + u[dof_z]
                )

                updated_coords[i] = Ferrite.Vec{3,Float64}(new_coords)
            end
            Ferrite.reinit!(facetvalues, cell, updated_coords, facet_idx)
            #fill!(fe_ext, 0.0) # Reset local force vector
            for qp in 1:getnquadpoints(facetvalues)
                dΓ = getdetJdV(facetvalues, qp)
                for i in 1:getnbasefunctions(facetvalues)
                    Nᵢ = shape_value(facetvalues, qp, i) # Shape function (vector for displacement)
                    # Add traction contribution to displacement block
                    fe_ext[BlockIndex((1), (i))] += (t_traction ⋅ Nᵢ) * dΓ
                end
            end
            # Assemble into global force vector using all DOFs (u and p)
            assemble!(F_ext, cell_dofs, fe_ext)
        end
    end
    return F_ext
end
####################################################################
function run_fem_hybrid(input::InputStruct)::RunResult
    cv_u = input.cell_values_u
    cv_p = input.cell_values_p
    fv = input.facet_values
    dh = input.dh
    ch = input.ch
    ΓN = input.facetsets
    tol = input.tol
    traction = input.tractions
    filename = input.filename
    output_dir = input.output_dir

    maxIterPerInc = input.maxIterPerInc
    totalTime = input.totalTime
    initInc = input.initInc
    minInc = input.minInc
    maxInc = input.maxInc
    totalInc = input.totalInc

    ndofs_dh = ndofs(dh)

    dof_U = input.dof_U
    dof_F = input.dof_F


    U = zeros(ndofs_dh)   # total displacement
    Uinit = zeros(ndofs_dh)   # running displacement for current step
    Residual = zeros(ndofs_dh)
    Total_F = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}()   # store displacements per step
    F_effect = Vector{Float64}()
    U_effect = Vector{Float64}()

    conv_flag = 0
    deltaT = initInc
    tot_time = 0.0
    tot_incr = 1
    failure_flag = 0

    while tot_time <= totalTime
        if tot_time == totalTime || deltaT == 0.0
            println("Analysis ended successfully.")
            break
        end
        if tot_incr > totalInc
            println("Maximum number of increments reached.")
            break
        end
        if failure_flag == 1
            break
        end

        n = totalTime / deltaT


        Incremental_F = zeros(ndofs_dh)

        # distribute traction incrementally
        traction_load = Dict{Int,Vector}(k => v ./ n for (k, v) in traction)
        F_ext = zeros(ndofs_dh)
        
        if input.model_type == :plane_stress
            assemble_traction_forces_twoD_hybrid!(F_ext, dh, ΓN, fv, traction_load, Uinit, cv_u, cv_p)
        elseif input.model_type == :plane_strain
            assemble_traction_forces_twoD_hybrid!(F_ext, dh, ΓN, fv, traction_load, Uinit, cv_u, cv_p)
        elseif input.model_type == :threeD
            assemble_traction_forces_threeD_hybrid!(F_ext, dh, ΓN, fv, traction_load, Uinit, cv_u, cv_p)
        else
            error("the model type is not correct")
        end

        Incremental_F = F_ext
        Total_F .+= Incremental_F


        for L = 1:maxIterPerInc

            K_nonlinear = allocate_matrix(dh)
            F_int = zeros(ndofs_dh)
            # assemble_global_plane_strain!(K_nonlinear, F_int, dh, cv, input, Uinit)
            assemble_global!(K_nonlinear, F_int, cv_u, cv_p, dh, input, Uinit)


            Residual .= Total_F - F_int
            Ferrite.update!(ch, tot_time + deltaT)  # Changed from tot_time
            apply!(K_nonlinear, Residual, ch)

            deltaU = K_nonlinear \ Residual

            Uinit .+= deltaU
            conv = norm(Residual) / (1 + norm(Total_F))
            if conv < tol
                conv_flag = 1
                break
            else
                conv_flag = 0
            end

        end

        if conv_flag == 1
            tot_time = tot_time + deltaT

            H1 = "Converged at increment $tot_incr, Δt = $deltaT, total time = $tot_time"
            println(H1)

            U = copy(Uinit)
            if isempty(dof_U) && isempty(dof_F)
                U_effect = []
                F_effect = []
            else
                U_mean = mean(U[dof_U])
                F_sum = sum(F_int[dof_F])
                push!(U_effect, U_mean)
                push!(F_effect, F_sum)
            end


            deltaT = deltaT * (1.25)  # increasing time increment for speed up
            if deltaT >= maxInc
                deltaT = maxInc
            end

            if deltaT >= (totalTime - tot_time)
                deltaT = (totalTime - tot_time)
            end

            tot_incr = tot_incr + 1

        elseif conv_flag == 0

            deltaT = deltaT / 4 # decreasing to improve convergence
            if deltaT < minInc
                println("Time integration has failed")
                failure_flag = 1
            end

        end
        push!(U_steps, copy(U))
    end
    VTKGridFile(joinpath(output_dir, filename), dh) do vtk
        write_solution(vtk, dh, U)
    end
    return RunResult(U_steps, U_effect, F_effect)
end