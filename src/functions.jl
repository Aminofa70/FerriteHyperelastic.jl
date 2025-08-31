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
"""
    initialize_solver([maxIterPerInc=500], [totalTime=1.0], [initInc=0.1],
                      [minInc=1e-5], [maxInc=0.2], [totalInc=500])

Initialize solver parameters for a time integration procedure.

# Arguments
- `maxIterPerInc::Int=500`: Maximum number of iterations allowed per increment.
- `totalTime::Float64=1.0`: Total simulation time.
- `initInc::Float64=0.1`: Initial increment size.
- `minInc::Float64=1e-5`: Minimum allowed increment size.
- `maxInc::Float64=0.2`: Maximum allowed increment size.
- `totalInc::Int=500`: Total number of increments.

# Behavior
- Ensures that `initInc` lies between `minInc` and `maxInc`.
- If `initInc > maxInc`, it is set to `maxInc`.
- If `initInc < minInc`, it is set to `minInc`.
- If `maxInc < minInc`, it is reset to `initInc`.
- Throws an error if the values are inconsistent after adjustment.

# Returns
A tuple:
(maxIterPerInc, totalTime, initInc, minInc, maxInc, totalInc)
"""
function initialize_solver(maxIterPerInc::Int=500, totalTime::Float64=1.0,
    initInc::Float64=0.1, minInc::Float64=1e-5,
    maxInc::Float64=0.2, totalInc::Int=500)

    # Ensure initInc is within bounds
    if initInc > maxInc
        initInc = maxInc
    end

    if initInc < minInc
        initInc = minInc
    end

    if maxInc < minInc
        maxInc = initInc
    end

    # Check consistency
    if initInc > maxInc || initInc < minInc
        error("Initial solver options are wrong.")
    end

    # Print info
    println("Solving the problem with the initial time integration of ",
        initInc, " and ", maxIterPerInc, " iterations per increment")

    return maxIterPerInc, totalTime, initInc, minInc, maxInc, totalInc
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
"""
function run_plane_strain(input::InputStruct)::RunResult
    cv         = input.cell_values
    fv         = input.facet_values
    dh         = input.dh
    ch         = input.ch
    ΓN         = input.facetsets
    tol        = input.tol
    traction   = input.tractions
    filename   = input.filename
    output_dir = input.output_dir
  
    maxIterPerInc = input.maxIterPerInc
    totalTime     = input.totalTime
    initInc       = input.initInc
    minInc        = input.minInc
    maxInc        = input.maxInc
    totalInc      = input.totalInc

    ndofs_dh = ndofs(dh)

    
    U      = zeros(ndofs_dh)   # total displacement
    Uinit  = zeros(ndofs_dh)   # running displacement for current step
    Residual = zeros(ndofs_dh)
    Total_F  = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}()   # store displacements per step


    conv_flag    = 0
    deltaT       = initInc 
    tot_time     = 0.0 
    tot_incr     = 1 
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
        traction_load = Dict{Int, Vector}(k => v ./ n for (k, v) in traction)
        F_ext = zeros(ndofs_dh)
        assemble_traction_forces_twoD!(F_ext, dh, ΓN, fv, traction_load, Uinit)
        
        Incremental_F = F_ext
        Total_F .+= Incremental_F
       
        
        for L = 1:maxIterPerInc
            
            K_nonlinear = allocate_matrix(dh)
            F_int = zeros(ndofs_dh)
            assemble_global_plane_strain!(K_nonlinear, F_int, dh, cv, input, Uinit)

            

            Residual .= Total_F - F_int
            Ferrite.update!(ch, tot_time + deltaT)  # Changed from tot_time
            apply!(K_nonlinear, Residual, ch)

            deltaU = K_nonlinear \ Residual
            
            Uinit .+= deltaU
            conv = norm(Residual) / (1 + norm(Total_F))
            if conv < tol
                conv_flag=1;
                break
            else
                conv_flag=0;
            end
            
        end 

        
        if conv_flag==1
            tot_time=tot_time+ deltaT;
    
            H1="CONVERGENCE IS REACHED FOR increment : $(tot_incr) with the converged time incremet of: $(deltaT) & total time of: $(tot_time)";println(H1);
            U=copy(Uinit);
    
            deltaT=deltaT*(1.25);  # increasing time increment for speed up
            if deltaT >= maxInc
                deltaT = maxInc;
            end
    
            if deltaT >= (totalTime - tot_time)
                deltaT=(totalTime - tot_time);
            end
    
            tot_incr=tot_incr+1;
    
        elseif conv_flag==0
            
            deltaT = deltaT/4; # decreasing to improve convergence
            if deltaT < minInc
                println("Time integration has failed")
                failure_flag=1;
            end
    
        end
        push!(U_steps, copy(U))
    end 
    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, U)
    end
    return RunResult(U_steps)
end
############################################################################################
############################################################################################
"""
"""
function run_plane_stress(input::InputStruct)::RunResult
    cv         = input.cell_values
    fv         = input.facet_values
    dh         = input.dh
    ch         = input.ch
    ΓN         = input.facetsets
    tol        = input.tol
    traction   = input.tractions
    filename   = input.filename
    output_dir = input.output_dir
  
    maxIterPerInc = input.maxIterPerInc
    totalTime     = input.totalTime
    initInc       = input.initInc
    minInc        = input.minInc
    maxInc        = input.maxInc
    totalInc      = input.totalInc

    ndofs_dh = ndofs(dh)

    
    U      = zeros(ndofs_dh)   # total displacement
    Uinit  = zeros(ndofs_dh)   # running displacement for current step
    Residual = zeros(ndofs_dh)
    Total_F  = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}()   # store displacements per step


    conv_flag    = 0
    deltaT       = initInc 
    tot_time     = 0.0 
    tot_incr     = 1 
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
        traction_load = Dict{Int, Vector}(k => v ./ n for (k, v) in traction)
        F_ext = zeros(ndofs_dh)
        assemble_traction_forces_twoD!(F_ext, dh, ΓN, fv, traction_load, Uinit)
        
        Incremental_F = F_ext
        Total_F .+= Incremental_F
       
        
        for L = 1:maxIterPerInc
            
            K_nonlinear = allocate_matrix(dh)
            F_int = zeros(ndofs_dh)
            assemble_global_plane_stress!(K_nonlinear, F_int, dh, cv, input, Uinit)

            

            Residual .= Total_F - F_int
            Ferrite.update!(ch, tot_time + deltaT)  # Changed from tot_time
            apply!(K_nonlinear, Residual, ch)

            deltaU = K_nonlinear \ Residual
            
            Uinit .+= deltaU
            conv = norm(Residual) / (1 + norm(Total_F))
            if conv < tol
                conv_flag=1;
                break
            else
                conv_flag=0;
            end
            
        end 

        
        if conv_flag==1
            tot_time=tot_time+ deltaT;
    
            H1="CONVERGENCE IS REACHED FOR increment : $(tot_incr) with the converged time incremet of: $(deltaT) & total time of: $(tot_time)";println(H1);
            U=copy(Uinit);
    
            deltaT=deltaT*(1.25);  # increasing time increment for speed up
            if deltaT >= maxInc
                deltaT = maxInc;
            end
    
            if deltaT >= (totalTime - tot_time)
                deltaT=(totalTime - tot_time);
            end
    
            tot_incr=tot_incr+1;
    
        elseif conv_flag==0
            
            deltaT = deltaT/4; # decreasing to improve convergence
            if deltaT < minInc
                println("Time integration has failed")
                failure_flag=1;
            end
    
        end
        push!(U_steps, copy(U))
    end 
    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, U)
    end
    return RunResult(U_steps)
end
############################################################################################
############################################################################################
"""
"""
function run_threeD(input::InputStruct)::RunResult
    cv         = input.cell_values
    fv         = input.facet_values
    dh         = input.dh
    ch         = input.ch
    ΓN         = input.facetsets
    tol        = input.tol
    traction   = input.tractions
    filename   = input.filename
    output_dir = input.output_dir
  
    maxIterPerInc = input.maxIterPerInc
    totalTime     = input.totalTime
    initInc       = input.initInc
    minInc        = input.minInc
    maxInc        = input.maxInc
    totalInc      = input.totalInc

    ndofs_dh = ndofs(dh)

    
    U      = zeros(ndofs_dh)   # total displacement
    Uinit  = zeros(ndofs_dh)   # running displacement for current step
    Residual = zeros(ndofs_dh)
    Total_F  = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}()   # store displacements per step


    conv_flag    = 0
    deltaT       = initInc 
    tot_time     = 0.0 
    tot_incr     = 1 
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
        traction_load = Dict{Int, Vector}(k => v ./ n for (k, v) in traction)
        F_ext = zeros(ndofs_dh)
        assemble_traction_forces_threeD!(F_ext, dh, ΓN, fv, traction_load, Uinit)
        
        Incremental_F = F_ext
        Total_F .+= Incremental_F
       
        
        for L = 1:maxIterPerInc
            
            K_nonlinear = allocate_matrix(dh)
            F_int = zeros(ndofs_dh)
            assemble_global_3D!(K_nonlinear, F_int, dh, cv, input, Uinit)

            

            Residual .= Total_F - F_int
            Ferrite.update!(ch, tot_time + deltaT)  # Changed from tot_time
            apply!(K_nonlinear, Residual, ch)

            deltaU = K_nonlinear \ Residual
            
            Uinit .+= deltaU
            conv = norm(Residual) / (1 + norm(Total_F))
            if conv < tol
                conv_flag=1;
                break
            else
                conv_flag=0;
            end
            
        end 

        
        if conv_flag==1
            tot_time=tot_time+ deltaT;
    
            H1="CONVERGENCE IS REACHED FOR increment : $(tot_incr) with the converged time incremet of: $(deltaT) & total time of: $(tot_time)";println(H1);
            U=copy(Uinit);
    
            deltaT=deltaT*(1.25);  # increasing time increment for speed up
            if deltaT >= maxInc
                deltaT = maxInc;
            end
    
            if deltaT >= (totalTime - tot_time)
                deltaT=(totalTime - tot_time);
            end
    
            tot_incr=tot_incr+1;
    
        elseif conv_flag==0
            
            deltaT = deltaT/4; # decreasing to improve convergence
            if deltaT < minInc
                println("Time integration has failed")
                failure_flag=1;
            end
    
        end
        push!(U_steps, copy(U))
    end 
    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, U)
    end
    return RunResult(U_steps)
end
############################################################################################
############################################################################################
"""

"""
function run_plane_strain_disp(input::InputStruct)::RunResult
    cv         = input.cell_values
    dh         = input.dh
    ch         = input.ch
    filename   = input.filename
    output_dir = input.output_dir
    
    tol          = input.tol
    displacement = input.displacement   # total prescribed displacement
    
    maxIterPerInc = input.maxIterPerInc
    totalTime     = input.totalTime
    initInc       = input.initInc
    minInc        = input.minInc
    maxInc        = input.maxInc
    totalInc      = input.totalInc

    ndofs_dh = ndofs(dh)
    U      = zeros(ndofs_dh)   # total displacement
    U_prev = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}()

    conv_flag    = 0
    deltaT       = initInc
    tot_time     = 0.0
    tot_incr     = 1
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

        # Adjust last increment to not overshoot
        if tot_time + deltaT > totalTime
            deltaT = totalTime - tot_time
        end

        # Target displacement factor
        λ = ((tot_time + deltaT) / totalTime) * displacement

        U_prev .= U
        conv_flag = 0

        # --- Newton loop ---
        for L = 1:maxIterPerInc
            K_nl = allocate_matrix(dh)
            F_int = zeros(ndofs_dh)
            assemble_global_plane_strain!(K_nl, F_int, dh, cv, input, U)

            Residual = -F_int

            Ferrite.update!(ch, λ)
            apply_zero!(K_nl, Residual, ch)

            ΔU = K_nl \ Residual
            U .+= ΔU
            apply!(U, ch)

            conv = norm(Residual[Ferrite.free_dofs(ch)]) / (1 + norm(F_int))
            if conv < tol
                conv_flag = 1
                break
            else
                conv_flag = 0
            end
        end
        # --- end Newton loop ---

        if conv_flag == 1
            tot_time = tot_time + deltaT

            println("CONVERGENCE IS REACHED FOR increment : $(tot_incr) with the converged time incremet of: $(deltaT) & total time of: $(tot_time)")

            push!(U_steps, copy(U))

            deltaT = deltaT * (1.25)  # increase step size
            if deltaT >= maxInc
                deltaT = maxInc
            end
            if deltaT >= (totalTime - tot_time)
                deltaT = (totalTime - tot_time)
            end

            tot_incr = tot_incr + 1

        elseif conv_flag == 0
            deltaT = deltaT / 4  # decrease step size
            U .= U_prev
            if deltaT < minInc
                println("Time integration has failed")
                failure_flag = 1
            end
        end
    end

    VTKGridFile(joinpath(output_dir, filename), dh) do vtk
        write_solution(vtk, dh, U)
    end
    return RunResult(U_steps)
end
############################################################################################
############################################################################################
"""

"""
function run_plane_stress_disp(input::InputStruct)::RunResult
    cv         = input.cell_values
    dh         = input.dh
    ch         = input.ch
    filename   = input.filename
    output_dir = input.output_dir
    
    tol          = input.tol
    displacement = input.displacement   # total prescribed displacement
    
    maxIterPerInc = input.maxIterPerInc
    totalTime     = input.totalTime
    initInc       = input.initInc
    minInc        = input.minInc
    maxInc        = input.maxInc
    totalInc      = input.totalInc

    ndofs_dh = ndofs(dh)
    U      = zeros(ndofs_dh)   # total displacement
    U_prev = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}()

    conv_flag    = 0
    deltaT       = initInc
    tot_time     = 0.0
    tot_incr     = 1
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

        # Adjust last increment to not overshoot
        if tot_time + deltaT > totalTime
            deltaT = totalTime - tot_time
        end

        # Target displacement factor
        λ = ((tot_time + deltaT) / totalTime) * displacement

        U_prev .= U
        conv_flag = 0

        # --- Newton loop ---
        for L = 1:maxIterPerInc
            K_nl = allocate_matrix(dh)
            F_int = zeros(ndofs_dh)
            assemble_global_plane_stress!(K_nl, F_int, dh, cv, input, U)

            Residual = -F_int

            Ferrite.update!(ch, λ)
            apply_zero!(K_nl, Residual, ch)

            ΔU = K_nl \ Residual
            U .+= ΔU
            apply!(U, ch)

            conv = norm(Residual[Ferrite.free_dofs(ch)]) / (1 + norm(F_int))
            if conv < tol
                conv_flag = 1
                break
            else
                conv_flag = 0
            end
        end
        # --- end Newton loop ---

        if conv_flag == 1
            tot_time = tot_time + deltaT

            println("CONVERGENCE IS REACHED FOR increment : $(tot_incr) with the converged time incremet of: $(deltaT) & total time of: $(tot_time)")

            push!(U_steps, copy(U))

            deltaT = deltaT * (1.25)   # increase step size
            if deltaT >= maxInc
                deltaT = maxInc
            end
            if deltaT >= (totalTime - tot_time)
                deltaT = (totalTime - tot_time)
            end

            tot_incr = tot_incr + 1

        elseif conv_flag == 0
            deltaT = deltaT / 4   # decrease step size
            U .= U_prev
            if deltaT < minInc
                println("Time integration has failed")
                failure_flag = 1
            end
        end
    end

    VTKGridFile(joinpath(output_dir, filename), dh) do vtk
        write_solution(vtk, dh, U)
    end
    return RunResult(U_steps)
end

############################################################################################
############################################################################################
"""
"""
function run_threeD_disp(input::InputStruct)::RunResult
    cv         = input.cell_values
    dh         = input.dh
    ch         = input.ch
    filename   = input.filename
    output_dir = input.output_dir
    
    tol          = input.tol
    displacement = input.displacement   # total prescribed displacement
    
    maxIterPerInc = input.maxIterPerInc
    totalTime     = input.totalTime
    initInc       = input.initInc
    minInc        = input.minInc
    maxInc        = input.maxInc
    totalInc      = input.totalInc

    ndofs_dh = ndofs(dh)
    U      = zeros(ndofs_dh)   # total displacement
    U_prev = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}()

    conv_flag    = 0
    deltaT       = initInc
    tot_time     = 0.0
    tot_incr     = 1
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

        # Adjust last increment to not overshoot
        if tot_time + deltaT > totalTime
            deltaT = totalTime - tot_time
        end

        # Target displacement factor
        λ = ((tot_time + deltaT) / totalTime) * displacement

        U_prev .= U
        conv_flag = 0

        # --- Newton loop ---
        for L = 1:maxIterPerInc
            K_nl = allocate_matrix(dh)
            F_int = zeros(ndofs_dh)
            assemble_global_3D!(K_nl, F_int, dh, cv, input, U)

            Residual = -F_int

            Ferrite.update!(ch, λ)
            apply_zero!(K_nl, Residual, ch)

            ΔU = K_nl \ Residual
            U .+= ΔU
            apply!(U, ch)

            conv = norm(Residual[Ferrite.free_dofs(ch)]) / (1 + norm(F_int))
            if conv < tol
                conv_flag = 1
                break
            else
                conv_flag = 0
            end
        end
        # --- end Newton loop ---

        if conv_flag == 1
            tot_time = tot_time + deltaT

            println("CONVERGENCE IS REACHED FOR increment : $(tot_incr) with the converged time incremet of: $(deltaT) & total time of: $(tot_time)")

            push!(U_steps, copy(U))

            deltaT = deltaT * (1.25)   # increase step size
            if deltaT >= maxInc
                deltaT = maxInc
            end
            if deltaT >= (totalTime - tot_time)
                deltaT = (totalTime - tot_time)
            end

            tot_incr = tot_incr + 1

        elseif conv_flag == 0
            deltaT = deltaT / 4   # decrease step size
            U .= U_prev
            if deltaT < minInc
                println("Time integration has failed")
                failure_flag = 1
            end
        end
    end

    VTKGridFile(joinpath(output_dir, filename), dh) do vtk
        write_solution(vtk, dh, U)
    end
    return RunResult(U_steps)
end

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