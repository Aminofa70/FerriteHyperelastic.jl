"""
    assemble_traction_forces_twoD!(
        F_ext,
        dh,
        facetsets::Vector,
        facetvalues,
        tractions::Dict{Int, <:AbstractVector},
        u::AbstractVector
    )

Assembles the global external force vector `F_ext` in 2D from prescribed
surface tractions applied on boundary facets.

Arguments
---------
- `F_ext` : Global external force vector (modified in place).
- `dh` : Discretization handler containing mesh topology, connectivity, and degrees of freedom.
- `facetsets::Vector` : List of facet sets (groups of boundary facets where tractions are applied).
- `facetvalues` : Precomputed quadrature and shape function values on facets.
- `tractions::Dict{Int, <:AbstractVector}` : Mapping from facet set ID to traction vectors
  (applied tractions in global coordinates).
- `u::AbstractVector` : Global displacement vector (may be used if tractions depend on current configuration).

Notes
-----
- For each facet in the provided `facetsets`, this function computes the
  traction contributions at quadrature points and assembles them into `F_ext`.
- The tractions can represent surface loads such as pressures or distributed forces
  acting on the model boundaries.
- This routine is intended for 2D problems; a corresponding 3D version would handle
  surface elements instead of line elements.
"""
function assemble_traction_forces_twoD!(
    F_ext,
    dh,
    facetsets::Vector,
    facetvalues,
    tractions::Dict{Int, <:AbstractVector},
    u::AbstractVector
)
    fe_ext = zeros(getnbasefunctions(facetvalues))
    for (idx, facetset) in enumerate(facetsets)
        t_traction = tractions[idx]
        for facetidx in facetset  # Iterate over FacetIndex objects in facetset
            # Get the cell and facet index
            cellid = facetidx[1]      # Cell ID
            facet_idx = facetidx[2]   # Local facet index
            cell = getcells(dh.grid)[cellid]
            
            # Get the original coordinates of the cell nodes
            cell_coords = getcoordinates(dh.grid, cellid)  # Vector{Vec{2, Float64}}, length 4
            
            # Update coordinates with displacement
            updated_coords = copy(cell_coords)
            cell_dofs = celldofs(dh, cellid)  # e.g., [ux1, uy1, ux2, uy2, ux3, uy3, ux4, uy4]
            for i in eachindex(cell_coords)# Iterate over nodes (1 to 4)
                dim = length(cell_coords[1])  # 2 for plane strain
                current_coords = cell_coords[i]
                # Get DOFs for the i-th node (x and y displacements)
                dof_x = cell_dofs[2*i-1]  # x-displacement DOF
                dof_y = cell_dofs[2*i]    # y-displacement DOF
                # Create new coordinates
                new_coords = (current_coords[1] + u[dof_x], current_coords[2] + u[dof_y])
                updated_coords[i] = Ferrite.Vec{2, Float64}(new_coords)
            end
            
            # Reinitialize FacetValues with updated coordinates
            reinit!(facetvalues, cell, updated_coords, facet_idx)
            
            # Compute traction forces
            fill!(fe_ext, 0.0)
            for qp in 1:getnquadpoints(facetvalues)
                dΓ = getdetJdV(facetvalues, qp)
                for i in 1:getnbasefunctions(facetvalues)
                    Nᵢ = shape_value(facetvalues, qp, i)
                    fe_ext[i] += t_traction ⋅ Nᵢ * dΓ
                end
            end
            assemble!(F_ext, cell_dofs, fe_ext)
        end
    end
    return F_ext
end
############################################################################################
############################################################################################
"""
    assemble_traction_forces_threeD!(
        F_ext,
        dh::DofHandler{3},
        facetsets::Vector,
        facetvalues::FacetValues,
        tractions::Dict{Int, <:AbstractVector},
        u::AbstractVector
    )

Assembles the global external force vector `F_ext` in 3D from prescribed
surface tractions applied on boundary facets.

Arguments
---------
- `F_ext` : Global external force vector (modified in place).
- `dh::DofHandler{3}` : Discretization handler for the 3D mesh, containing topology,
  connectivity, and degrees of freedom.
- `facetsets::Vector` : List of facet sets (groups of boundary surface facets where tractions are applied).
- `facetvalues::FacetValues` : Quadrature and shape function values on the facets.
- `tractions::Dict{Int, <:AbstractVector}` : Mapping from facet set ID to traction vectors
  (applied tractions in global 3D coordinates).
- `u::AbstractVector` : Global displacement vector (may be used if tractions depend on the current configuration).

Notes
-----
- For each facet in the provided `facetsets`, this function computes the
  traction contributions at quadrature points on the surface and assembles them into `F_ext`.
- Typical use cases include pressures, surface loads, or other distributed forces
  acting on the boundary of a 3D domain.
- This routine is the 3D analogue of `assemble_traction_forces_twoD!`.
"""
function assemble_traction_forces_threeD!(
    F_ext,
    dh::DofHandler{3},
    facetsets::Vector,
    facetvalues::FacetValues,
    tractions::Dict{Int, <:AbstractVector},
    u::AbstractVector
)
    fe_ext = zeros(getnbasefunctions(facetvalues))
    for (idx, facetset) in enumerate(facetsets)
        t_traction = tractions[idx]  # Vec{3, Float64} for 3D
        for facetidx in facetset  # Iterate over FacetIndex objects in facetset
            # Get the cell and facet index
            cellid = facetidx[1]      # Cell ID
            facet_idx = facetidx[2]   # Local facet index
            cell = getcells(dh.grid)[cellid]
            
            # Get the original coordinates of the cell nodes
            cell_coords = getcoordinates(dh.grid, cellid)  # Vector{Vec{3, Float64}}
            
            # Update coordinates with displacement
            updated_coords = copy(cell_coords)
            cell_dofs = celldofs(dh, cellid)  # e.g., [ux1, uy1, uz1, ux2, uy2, uz2, ...]
            for i in eachindex(cell_coords) # Iterate over nodes (e.g., 8 for hexahedron)
                dim = length(cell_coords[1])  # 3 for 3D
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
                updated_coords[i] = Ferrite.Vec{3, Float64}(new_coords)
            end
            
            # Reinitialize FacetValues with updated coordinates
            reinit!(facetvalues, cell, updated_coords, facet_idx)
            
            # Compute traction forces
            fill!(fe_ext, 0.0)
            for qp in 1:getnquadpoints(facetvalues)
                dΓ = getdetJdV(facetvalues, qp)
                for i in 1:getnbasefunctions(facetvalues)
                    Nᵢ = shape_value(facetvalues, qp, i)  # Scalar shape function
                    fe_ext[i] += t_traction ⋅ Nᵢ * dΓ    # Traction is Vec{3, Float64}
                end
            end
            assemble!(F_ext, cell_dofs, fe_ext)
        end
    end
    return F_ext
end
############################################################################################
############################################################################################
struct RunResult
    U_steps::Vector{Vector{Float64}}  # displacement vectors at each converged step
end
####
############################################################################################
############################################################################################
"""
    run_plane_strain(input::InputStruct) :: RunResult

Runs a finite element simulation under plane strain conditions and returns
the displacement field at each load/time step.

Arguments
---------
- `input::InputStruct` : User-defined input data containing mesh information,
  boundary conditions, material properties, solver settings, and load steps.

Returns
-------
- `RunResult` : Object containing the displacement vector(s) for all
  increments/steps of the simulation.

Notes
-----
- The function sets up the mesh, degrees of freedom, and solver based on `input`.
- For each step, the global system is assembled via
  `assemble_global_plane_strain!`, external forces are applied,
  and the displacement solution is computed.
- Only displacements are stored in `RunResult`; stresses, strains,
  or reactions must be post-processed separately.
- Plane strain assumption is enforced: ε₃₃ = 0 but σ₃₃ ≠ 0.
"""
function run_plane_strain(input::InputStruct)::RunResult
    cv        = input.cell_values
    fv        = input.facet_values
    dh        = input.dh
    ch        = input.ch
    ΓN        = input.facetsets
    tol       = input.tol
    n_load_steps = input.n_load_steps
    traction  = input.tractions
    filename  = input.filename
    output_dir= input.output_dir
    n_iter_NR = input.n_iter_NR

    ndofs_dh = ndofs(dh)

    U      = zeros(ndofs_dh)   # total displacement
    Uinit  = zeros(ndofs_dh)   # running displacement for current step
    Residual = zeros(ndofs_dh)
    Total_F  = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}()   # store displacements per step

    for O in 1:n_load_steps
        Incremental_F = zeros(ndofs_dh)

        # distribute traction incrementally
        traction_load = Dict{Int, Vector}(k => v ./ n_load_steps for (k, v) in traction)
        F_ext = zeros(ndofs_dh)
        assemble_traction_forces_twoD!(F_ext, dh, ΓN, fv, traction_load, Uinit)

        Incremental_F .+= F_ext
        Total_F .+= Incremental_F

        conv_flag = false
        conv = Inf

        # --- Newton iterations ---
        for L in 1:n_iter_NR
            K_nonlinear = allocate_matrix(dh)
            F_int = zeros(ndofs_dh)
            assemble_global_plane_strain!(K_nonlinear, F_int, dh, cv, input, Uinit)

            Residual .= Total_F - F_int
            Ferrite.update!(ch, O)
            apply!(K_nonlinear, Residual, ch)

            deltaU = K_nonlinear \ Residual
            Uinit .+= deltaU

            conv = norm(Residual) / (1 + norm(Total_F))
            if conv < tol
                conv_flag = true
                break
            end
        end

        if !conv_flag
            println("❌ Convergence NOT reached at step $O")
            break
        else
            println("✅ Convergence reached at step $O")
            U .= Uinit
        end

        

        # Store displacement vector
        push!(U_steps, copy(U))
    end
    # Write VTK for this step
    
    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, U)
    end

    return RunResult(U_steps)
end
############################################################################################
############################################################################################
"""
    run_plane_stress(input::InputStruct) :: RunResult

Runs a finite element simulation under plane stress conditions and returns
the displacement field at each load/time step.

Arguments
---------
- `input::InputStruct` : User-defined input data containing mesh information,
  boundary conditions, material properties, solver settings, and load steps.

Returns
-------
- `RunResult` : Object containing the displacement vector(s) for all
  increments/steps of the simulation.

Notes
-----
- The function sets up the mesh, degrees of freedom, and solver based on `input`.
- For each step, the global system is assembled via
  `assemble_global_plane_stress!`, external forces are applied,
  and the displacement solution is computed.
- Only displacements are stored in `RunResult`; stresses, strains,
  or reactions must be post-processed separately.
- Plane stress condition is enforced: σ₃₃ = 0, with the out-of-plane stretch λ₃
  solved iteratively using `solve_lambda3`.
"""
function run_plane_stress(input::InputStruct)::RunResult
    cv        = input.cell_values
    fv        = input.facet_values
    dh        = input.dh
    ch        = input.ch
    ΓN        = input.facetsets
    tol       = input.tol
    n_load_steps = input.n_load_steps
    traction  = input.tractions
    filename  = input.filename
    output_dir= input.output_dir
    n_iter_NR = input.n_iter_NR

    ndofs_dh = ndofs(dh)

    U      = zeros(ndofs_dh)   # total displacement
    Uinit  = zeros(ndofs_dh)   # running displacement for current step
    Residual = zeros(ndofs_dh)
    Total_F  = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}()   # store displacements per step

    for O in 1:n_load_steps
        Incremental_F = zeros(ndofs_dh)

        # distribute traction incrementally
        traction_load = Dict{Int, Vector}(k => v ./ n_load_steps for (k, v) in traction)
        F_ext = zeros(ndofs_dh)
        assemble_traction_forces_twoD!(F_ext, dh, ΓN, fv, traction_load, Uinit)

        Incremental_F .+= F_ext
        Total_F .+= Incremental_F

        conv_flag = false
        conv = Inf

        # --- Newton iterations ---
        for L in 1:n_iter_NR
            K_nonlinear = allocate_matrix(dh)
            F_int = zeros(ndofs_dh)
            assemble_global_plane_stress!(K_nonlinear, F_int, dh, cv, input, Uinit)

            Residual .= Total_F - F_int
            Ferrite.update!(ch, O)
            apply!(K_nonlinear, Residual, ch)

            deltaU = K_nonlinear \ Residual
            Uinit .+= deltaU

            conv = norm(Residual) / (1 + norm(Total_F))
            if conv < tol
                conv_flag = true
                break
            end
        end

        if !conv_flag
            println("❌ Convergence NOT reached at step $O")
            break
        else
            println("✅ Convergence reached at step $O")
            U .= Uinit
        end

        

        # Store displacement vector
        push!(U_steps, copy(U))
    end
    # Write VTK for this step
    
    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, U)
    end

    return RunResult(U_steps)
end
############################################################################################
############################################################################################
"""
    run_threeD(input::InputStruct) :: RunResult

Runs a finite element simulation in full 3D and returns
the displacement field at each load/time step.

Arguments
---------
- `input::InputStruct` : User-defined input data containing mesh information,
  boundary conditions, material properties, solver settings, and load steps.

Returns
-------
- RunResult` : Object containing the displacement vector(s) for all
  increments/steps of the simulation.`

Notes
-----
- The function sets up the mesh, degrees of freedom, and solver based on `input`.
- For each step, the global system is assembled via
  `assemble_global_3D!`, external forces are applied,
  and the displacement solution is computed.
- Only displacements are stored in `RunResult`; stresses, strains,
  or reactions must be post-processed separately.
- Full 3D formulation is used: all three displacement components
  and stress/strain components are active (no plane stress/strain assumptions).
"""
function run_threeD(input::InputStruct)::RunResult
    cv        = input.cell_values
    fv        = input.facet_values
    dh        = input.dh
    ch        = input.ch
    ΓN        = input.facetsets
    tol       = input.tol
    n_load_steps = input.n_load_steps
    traction  = input.tractions
    filename  = input.filename
    output_dir= input.output_dir
    n_iter_NR = input.n_iter_NR

    ndofs_dh = ndofs(dh)

    U      = zeros(ndofs_dh)   # total displacement
    Uinit  = zeros(ndofs_dh)   # running displacement for current step
    Residual = zeros(ndofs_dh)
    Total_F  = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}()   # store displacements per step

    for O in 1:n_load_steps
        Incremental_F = zeros(ndofs_dh)

        # distribute traction incrementally
        traction_load = Dict{Int, Vector}(k => v ./ n_load_steps for (k, v) in traction)
        F_ext = zeros(ndofs_dh)
        assemble_traction_forces_threeD!(F_ext, dh, ΓN, fv, traction_load, Uinit)

        Incremental_F .+= F_ext
        Total_F .+= Incremental_F

        conv_flag = false
        conv = Inf

        # --- Newton iterations ---
        for L in 1:n_iter_NR
            K_nonlinear = allocate_matrix(dh)
            F_int = zeros(ndofs_dh)
            assemble_global_3D!(K_nonlinear, F_int, dh, cv, input, Uinit)

            Residual .= Total_F - F_int
            Ferrite.update!(ch, O)
            apply!(K_nonlinear, Residual, ch)

            deltaU = K_nonlinear \ Residual
            Uinit .+= deltaU

            conv = norm(Residual) / (1 + norm(Total_F))
            if conv < tol
                conv_flag = true
                break
            end
        end

        if !conv_flag
            println("❌ Convergence NOT reached at step $O")
            break
        else
            println("✅ Convergence reached at step $O")
            U .= Uinit
        end

        

        # Store displacement vector
        push!(U_steps, copy(U))
    end
    # Write VTK for this step
    
    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, U)
    end

    return RunResult(U_steps)
end
############################################################################################
############################################################################################
"""
    run_plane_strain_disp(input::InputStruct) -> RunResult

Solve a nonlinear plane strain finite element problem using incremental
displacement loading and a Newton–Raphson scheme with line search. 

The solver:
  • Initializes global vectors and stiffness matrix.  
  • Applies displacement increments over `n_load_steps`.  
  • At each step, assembles residual and tangent stiffness, solves for
    displacement updates, and checks convergence against `tol`.  
  • Uses adaptive load stepping if convergence fails and line search to
    improve robustness.  
  • Writes the final displacement field to a VTK file for visualization.  

Returns a `RunResult` containing the converged displacement solution.
"""
function run_plane_strain_disp(input::InputStruct)::RunResult
    cv        = input.cell_values
    dh        = input.dh
    ch        = input.ch
    filename  = input.filename
    output_dir= input.output_dir
    n_load_steps = input.n_load_steps
    tol       = input.tol
    displacement  = input.displacement
    n_iter_NR = input.n_iter_NR
    LINE_SEARCH_STEPS = input.LINE_SEARCH_STEPS
    # Get number of degrees of freedom
    _ndofs = ndofs(dh)
    
    # Pre-allocation of vectors for the solution and Newton increments
    w = zeros(_ndofs)      # Solution vector (displacements and possibly pressure)
    ΔΔw = zeros(_ndofs)    # Newton increment vector
    w_prev = zeros(_ndofs) # Store previous solution for adaptive time stepping
    
    # Create the sparse matrix and residual vector
    K = allocate_matrix(dh)
    f = zeros(_ndofs)


    Δλ = displacement/ n_load_steps  # Incremental load factor
    
    U_steps = Vector{Vector{Float64}}()   # store displacements per step

    # Incremental displacement loop
    for step in 1:n_load_steps
        # Current load factor (scales displacement boundary conditions)
        λ = step * Δλ
        
        # Update constraint handler for current load factor
        Ferrite.update!(ch, λ)
        
        # Initialize solution with scaled Dirichlet boundary conditions
        apply!(w, ch)
        
        newton_itr = 0
        norm_res = Inf
        #prog = ProgressMeter.ProgressThresh(NEWTON_TOL; desc = "Solving @ step $step of $NUM_STEPS (λ = $λ);")
        fill!(ΔΔw, 0.0)
        w_prev .= w  # Store current solution

        # Newton iteration loop with line search
        while true
            newton_itr += 1
            
            # Assemble stiffness matrix and internal force vector
            K_nonlinear = allocate_matrix(dh)
            F_int = zeros(_ndofs)
            assemble_global_plane_strain!(K_nonlinear, F_int, dh, cv, input, w)

            # Compute residual
            K .= K_nonlinear
            f .= F_int
            
            # Apply Dirichlet boundary conditions to residual and stiffness matrix
            apply_zero!(K, f, ch)
            
            # Compute residual norm for free degrees of freedom
            norm_res = norm(f[Ferrite.free_dofs(ch)])
            
            # # Update progress meter
            # ProgressMeter.update!(prog, norm_res; showvalues = [(:iter, newton_itr)])
            
            # Check for convergence
            if norm_res < tol
                #@info "Converged at step $step (λ = $λ) after $newton_itr iterations"
                println("✅ Convergence reached at step $step")
                break
            elseif newton_itr > n_iter_NR
                @warn "Reached maximum Newton iterations at step $step (λ = $λ), attempting smaller increment"
                
                # Reduce load factor increment and retry
                Δλ /= 2
                if Δλ < 1e-6
                    error("Load factor increment too small, aborting at step $step")
                end
                step -= 1  # Step back
                λ = step * Δλ
                w .= w_prev  # Restore previous solution
                Ferrite.update!(ch, λ)
                apply!(w, ch)
                newton_itr = 0
                fill!(ΔΔw, 0.0)
                continue
            end
            
            # Solve for Newton increment
            ΔΔw .= K \ f
            apply_zero!(ΔΔw, ch)
            
            # Simple line search
            α = 1.0  # Line search parameter
            w_temp = copy(w)
            for ls in 1:LINE_SEARCH_STEPS
                w_temp .= w .- α .* ΔΔw
                F_int_temp = zeros(_ndofs)
                assemble_global_plane_strain!(K_nonlinear, F_int_temp, dh, cv, input, w_temp)
                apply_zero!(K_nonlinear, F_int_temp, ch)
                norm_res_temp = norm(F_int_temp[Ferrite.free_dofs(ch)])
                
                if norm_res_temp < norm_res
                    w .= w_temp
                    break
                end
                α /= 2
                if ls == LINE_SEARCH_STEPS
                    w .= w .- ΔΔw
                end
            end
        end
        
        # Apply Dirichlet boundary conditions to solution after convergence
        apply!(w, ch)
        
        push!(U_steps, copy(w))
      
    end

    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, w)
    end
    return RunResult(U_steps)
end

############################################################################################
############################################################################################
"""
    run_plane_stress_disp(input::InputStruct) -> RunResult

Solve a nonlinear plane stress finite element problem using incremental
displacement loading and a Newton–Raphson scheme with line search. 

The solver:
  • Initializes global vectors and stiffness matrix.  
  • Applies displacement increments over `n_load_steps`.  
  • At each step, assembles residual and tangent stiffness, solves for
    displacement updates, and checks convergence against `tol`.  
  • Uses adaptive load stepping if convergence fails and line search to
    improve robustness.  
  • Writes the final displacement field to a VTK file for visualization.  

Returns a `RunResult` containing the converged displacement solution.
"""
function run_plane_stress_disp(input::InputStruct)::RunResult
    cv        = input.cell_values
    dh        = input.dh
    ch        = input.ch
    filename  = input.filename
    output_dir= input.output_dir
    n_load_steps = input.n_load_steps
    tol       = input.tol
    displacement  = input.displacement
    n_iter_NR = input.n_iter_NR
    LINE_SEARCH_STEPS = input.LINE_SEARCH_STEPS
    # Get number of degrees of freedom
    _ndofs = ndofs(dh)
    
    # Pre-allocation of vectors for the solution and Newton increments
    w = zeros(_ndofs)      # Solution vector (displacements and possibly pressure)
    ΔΔw = zeros(_ndofs)    # Newton increment vector
    w_prev = zeros(_ndofs) # Store previous solution for adaptive time stepping
    
    # Create the sparse matrix and residual vector
    K = allocate_matrix(dh)
    f = zeros(_ndofs)


    Δλ = displacement/ n_load_steps  # Incremental load factor
    
    U_steps = Vector{Vector{Float64}}()   # store displacements per step

    # Incremental displacement loop
    for step in 1:n_load_steps
        # Current load factor (scales displacement boundary conditions)
        λ = step * Δλ
        
        # Update constraint handler for current load factor
        Ferrite.update!(ch, λ)
        
        # Initialize solution with scaled Dirichlet boundary conditions
        apply!(w, ch)
        
        newton_itr = 0
        norm_res = Inf
        #prog = ProgressMeter.ProgressThresh(NEWTON_TOL; desc = "Solving @ step $step of $NUM_STEPS (λ = $λ);")
        fill!(ΔΔw, 0.0)
        w_prev .= w  # Store current solution

        # Newton iteration loop with line search
        while true
            newton_itr += 1
            
            # Assemble stiffness matrix and internal force vector
            K_nonlinear = allocate_matrix(dh)
            F_int = zeros(_ndofs)
            assemble_global_plane_stress!(K_nonlinear, F_int, dh, cv, input, w)

            # Compute residual
            K .= K_nonlinear
            f .= F_int
            
            # Apply Dirichlet boundary conditions to residual and stiffness matrix
            apply_zero!(K, f, ch)
            
            # Compute residual norm for free degrees of freedom
            norm_res = norm(f[Ferrite.free_dofs(ch)])
            
            # # Update progress meter
            # ProgressMeter.update!(prog, norm_res; showvalues = [(:iter, newton_itr)])
            
            # Check for convergence
            if norm_res < tol
                #@info "Converged at step $step (λ = $λ) after $newton_itr iterations"
                println("✅ Convergence reached at step $step")
                break
            elseif newton_itr > n_iter_NR
                @warn "Reached maximum Newton iterations at step $step (λ = $λ), attempting smaller increment"
                
                # Reduce load factor increment and retry
                Δλ /= 2
                if Δλ < 1e-6
                    error("Load factor increment too small, aborting at step $step")
                end
                step -= 1  # Step back
                λ = step * Δλ
                w .= w_prev  # Restore previous solution
                Ferrite.update!(ch, λ)
                apply!(w, ch)
                newton_itr = 0
                fill!(ΔΔw, 0.0)
                continue
            end
            
            # Solve for Newton increment
            ΔΔw .= K \ f
            apply_zero!(ΔΔw, ch)
            
            # Simple line search
            α = 1.0  # Line search parameter
            w_temp = copy(w)
            for ls in 1:LINE_SEARCH_STEPS
                w_temp .= w .- α .* ΔΔw
                F_int_temp = zeros(_ndofs)
                assemble_global_plane_stress!(K_nonlinear, F_int_temp, dh, cv, input, w_temp)
                apply_zero!(K_nonlinear, F_int_temp, ch)
                norm_res_temp = norm(F_int_temp[Ferrite.free_dofs(ch)])
                
                if norm_res_temp < norm_res
                    w .= w_temp
                    break
                end
                α /= 2
                if ls == LINE_SEARCH_STEPS
                    w .= w .- ΔΔw
                end
            end
        end
        
        # Apply Dirichlet boundary conditions to solution after convergence
        apply!(w, ch)
        
        push!(U_steps, copy(w))
      
    end

    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, w)
    end
    return RunResult(U_steps)
end

############################################################################################
############################################################################################

function run_threeD_disp(input::InputStruct)::RunResult
    cv        = input.cell_values
    dh        = input.dh
    ch        = input.ch
    filename  = input.filename
    output_dir= input.output_dir
    n_load_steps = input.n_load_steps
    tol       = input.tol
    displacement  = input.displacement
    n_iter_NR = input.n_iter_NR
    LINE_SEARCH_STEPS = input.LINE_SEARCH_STEPS
    # Get number of degrees of freedom
    _ndofs = ndofs(dh)
    
    # Pre-allocation of vectors for the solution and Newton increments
    w = zeros(_ndofs)      # Solution vector (displacements and possibly pressure)
    ΔΔw = zeros(_ndofs)    # Newton increment vector
    w_prev = zeros(_ndofs) # Store previous solution for adaptive time stepping
    
    # Create the sparse matrix and residual vector
    K = allocate_matrix(dh)
    f = zeros(_ndofs)


    Δλ = displacement/ n_load_steps  # Incremental load factor
    
    U_steps = Vector{Vector{Float64}}()   # store displacements per step

    # Incremental displacement loop
    for step in 1:n_load_steps
        # Current load factor (scales displacement boundary conditions)
        λ = step * Δλ
        
        # Update constraint handler for current load factor
        Ferrite.update!(ch, λ)
        
        # Initialize solution with scaled Dirichlet boundary conditions
        apply!(w, ch)
        
        newton_itr = 0
        norm_res = Inf
        #prog = ProgressMeter.ProgressThresh(NEWTON_TOL; desc = "Solving @ step $step of $NUM_STEPS (λ = $λ);")
        fill!(ΔΔw, 0.0)
        w_prev .= w  # Store current solution

        # Newton iteration loop with line search
        while true
            newton_itr += 1
            
            # Assemble stiffness matrix and internal force vector
            K_nonlinear = allocate_matrix(dh)
            F_int = zeros(_ndofs)
            assemble_global_3D!(K_nonlinear, F_int, dh, cv, input, w)

            # Compute residual
            K .= K_nonlinear
            f .= F_int
            
            # Apply Dirichlet boundary conditions to residual and stiffness matrix
            apply_zero!(K, f, ch)
            
            # Compute residual norm for free degrees of freedom
            norm_res = norm(f[Ferrite.free_dofs(ch)])
            
            # # Update progress meter
            # ProgressMeter.update!(prog, norm_res; showvalues = [(:iter, newton_itr)])
            
            # Check for convergence
            if norm_res < tol
                #@info "Converged at step $step (λ = $λ) after $newton_itr iterations"
                println("✅ Convergence reached at step $step")
                break
            elseif newton_itr > n_iter_NR
                @warn "Reached maximum Newton iterations at step $step (λ = $λ), attempting smaller increment"
                
                # Reduce load factor increment and retry
                Δλ /= 2
                if Δλ < 1e-6
                    error("Load factor increment too small, aborting at step $step")
                end
                step -= 1  # Step back
                λ = step * Δλ
                w .= w_prev  # Restore previous solution
                Ferrite.update!(ch, λ)
                apply!(w, ch)
                newton_itr = 0
                fill!(ΔΔw, 0.0)
                continue
            end
            
            # Solve for Newton increment
            ΔΔw .= K \ f
            apply_zero!(ΔΔw, ch)
            
            # Simple line search
            α = 1.0  # Line search parameter
            w_temp = copy(w)
            for ls in 1:LINE_SEARCH_STEPS
                w_temp .= w .- α .* ΔΔw
                F_int_temp = zeros(_ndofs)
                assemble_global_3D!(K_nonlinear, F_int_temp, dh, cv, input, w_temp)
                apply_zero!(K_nonlinear, F_int_temp, ch)
                norm_res_temp = norm(F_int_temp[Ferrite.free_dofs(ch)])
                
                if norm_res_temp < norm_res
                    w .= w_temp
                    break
                end
                α /= 2
                if ls == LINE_SEARCH_STEPS
                    w .= w .- ΔΔw
                end
            end
        end
        
        # Apply Dirichlet boundary conditions to solution after convergence
        apply!(w, ch)
        
        push!(U_steps, copy(w))
      
    end

    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, w)
    end
    return RunResult(U_steps)
end

############################################################################################
############################################################################################
function run_fem(input::InputStruct)
    if input.load_type == :traction
        if input.model_type == :plane_stress
            run_plane_stress(input)
        elseif input.model_type == :plane_strain
            run_plane_strain(input)
        elseif input.model_type == :threeD
            run_threeD(input)
        else
            error("Unknown model_type: $(input.model_type)")
        end

    elseif input.load_type == :displacement
        if input.model_type == :plane_stress
            run_plane_stress_disp(input)
        elseif input.model_type == :plane_strain
            run_plane_strain_disp(input)
        elseif input.model_type == :threeD
            run_threeD_disp(input)
        else
            error("Unknown model_type: $(input.model_type)")
        end

    else
        error("Unknown load_type: $(input.load_type)")
    end
end

############################################################################################
############################################################################################