############################################################################################
############################################################################################
"""
    initialize_solver([maxIterPerInc=500], [totalTime=1.0], [initInc=0.1],
                      [minInc=1e-5], [maxInc=0.2], [totalInc=500])

Initialize solver parameters for a time integration procedure.
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
    U_steps::Vector{Vector{Float64}} 
    U_effect::Vector{Float64} 
    F_effect::Vector{Float64} 
end

####
############################################################################################
############################################################################################
"""
    run_plane_strain(input)

run finite element for 2D plane strain
"""
function run_plane_strain(input::InputStruct)
    # run_plane_strain(input::InputStruct)::RunResult
    cv         = input.cell_values
    fv         = input.facet_values
    dh         = input.dh
    ch         = input.ch
    ΓN         = input.facetsets
    tol        = input.tol
    traction   = input.tractions
    filename   = input.filename
    output_dir = input.output_dir

    dof_U = input.dof_U
    dof_F = input.dof_F
  
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

    

    
    F_int = zeros(ndofs_dh)  
     
    U_steps = Vector{Vector{Float64}}() 
    F_effect = Vector{Float64}()   
    U_effect = Vector{Float64}() 
      
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
        
        #assemble_traction_forces_twoD!(F_ext, dh, ΓN, fv, traction_load, Uinit)
        assemble_traction_forces_plane_strain_twoD!(F_ext, dh, ΓN, fv, traction_load, Uinit)

        
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

            apply_zero!(deltaU, ch)
            
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

            tot_time =  tot_time+ deltaT;

            # H1 = "CONVERGENCE IS REACHED FOR increment: $tot_incr " *
            #      "with the converged time increment of: $deltaT " *
            #      "& total time of: $tot_time"

            H1 = "Converged at increment $tot_incr, Δt = $deltaT, total time = $tot_time"
            println(H1)

            println(H1);
            U = copy(Uinit);
            # find the mean of U

            if isempty(dof_U) && isempty(dof_F)
                U_effect = []
                F_effect = []
            else
                U_mean = mean(U[dof_U])
                F_sum = sum(F_int[dof_F])
                push!(U_effect, U_mean)
                push!(F_effect, F_sum)
            end
    
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
        # push!(U_effect, U_mean )
        # push!(F_effect, F_sum)
    end 
    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, U)
    end
  
    
    return RunResult(U_steps, U_effect, F_effect)
end
############################################################################################
############################################################################################
"""
    run_plane_stress(input)
   
run finite element for plane stress
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

    dof_U = input.dof_U
    dof_F = input.dof_F

    U      = zeros(ndofs_dh)   # total displacement
    Uinit  = zeros(ndofs_dh)   # running displacement for current step
    Residual = zeros(ndofs_dh)
    Total_F  = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}() 
    F_effect = Vector{Float64}()   
    U_effect = Vector{Float64}() 


    conv_flag    = 0
    deltaT       = initInc 
    tot_time     = 0.0 
    tot_incr     = 1 
    failure_flag = 0 
    
    while tot_time <= totalTime
        
        if tot_time == totalTime || deltaT < 1e-15
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
        assemble_traction_forces_plane_stress_twoD!(F_ext, dh, ΓN, fv, traction_load, Uinit, input)
        
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
            apply_zero!(deltaU, ch)

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
    
            H1 = "Converged at increment $tot_incr, Δt = $deltaT, total time = $tot_time"
            println(H1);
            U=copy(Uinit);

            if isempty(dof_U) && isempty(dof_F)
                U_effect = []
                F_effect = []
            else
                U_mean = mean(U[dof_U])
                F_sum = sum(F_int[dof_F])
                push!(U_effect, U_mean)
                push!(F_effect, F_sum)
            end
    
            deltaT=deltaT*(1.25);  # increasing time increment for speed up
            if deltaT >= maxInc
                deltaT = maxInc;
            end
    
            if deltaT >= (totalTime - tot_time)
                deltaT=0.0;
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
        # push!(U_effect, U_mean )
        # push!(F_effect, F_sum)
    end 
    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, U)
    end
    return RunResult(U_steps, U_effect, F_effect)
end
############################################################################################
############################################################################################
"""
    run_threeD(input)


run finite element for 3D 
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

    dof_U = input.dof_U
    dof_F = input.dof_F

    U      = zeros(ndofs_dh)   # total displacement
    Uinit  = zeros(ndofs_dh)   # running displacement for current step
    Residual = zeros(ndofs_dh)
    Total_F  = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}() 
    F_effect = Vector{Float64}()   
    U_effect = Vector{Float64}() 


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
            
            apply_zero!(deltaU, ch)

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
    
            H1 = "Converged at increment $tot_incr, Δt = $deltaT, total time = $tot_time"
            println(H1);
            U=copy(Uinit);
            if isempty(dof_U) && isempty(dof_F)
                U_effect = []
                F_effect = []
            else
                U_mean = mean(U[dof_U])
                F_sum = sum(F_int[dof_F])
                push!(U_effect, U_mean)
                push!(F_effect, F_sum)
            end

            deltaT = deltaT*(1.25);  # increasing time increment for speed up
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
        # push!(U_effect, U_mean )
        # push!(F_effect, F_sum)
    end 
    VTKGridFile(joinpath(output_dir,filename), dh) do vtk
        write_solution(vtk, dh, U)
    end
    return RunResult(U_steps, U_effect, F_effect)
end
############################################################################################
############################################################################################
"""
     run_plane_strain_disp(input)

run finite element for plane stain with displacement load
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

    dof_U = input.dof_U
    dof_F = input.dof_F


    ndofs_dh = ndofs(dh)
    U      = zeros(ndofs_dh)   # total displacement
    U_prev = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}() 
    F_effect = Vector{Float64}()   
    U_effect = Vector{Float64}() 


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
            apply_zero!(ΔU, ch)
            U .+= ΔU
            #apply!(U, ch)

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

            H1 = "Converged at increment $tot_incr, Δt = $deltaT, total time = $tot_time"
            println(H1);
            if isempty(dof_U) && isempty(dof_F)
                U_effect = []
                F_effect = []
            else
                U_mean = mean(U[dof_U])
                F_sum = sum(F_int[dof_F])
                push!(U_effect, U_mean)
                push!(F_effect, F_sum)
            end

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
    return RunResult(U_steps, U_effect, F_effect)
end
############################################################################################
############################################################################################
"""
    run_plane_stress_disp(input)

run finite element for plane stress with displacement load 
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

    dof_U = input.dof_U
    dof_F = input.dof_F

    U      = zeros(ndofs_dh)   # total displacement
    U_prev = zeros(ndofs_dh)

    U_steps = Vector{Vector{Float64}}() 
    F_effect = Vector{Float64}()   
    U_effect = Vector{Float64}() 

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
            apply_zero!(ΔU, ch)
            U .+= ΔU
            #apply!(U, ch)

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

            H1 = "Converged at increment $tot_incr, Δt = $deltaT, total time = $tot_time"
            println(H1)

            if isempty(dof_U) && isempty(dof_F)
                U_effect = []
                F_effect = []
            else
                U_mean = mean(U[dof_U])
                F_sum = sum(F_int[dof_F])
                push!(U_effect, U_mean)
                push!(F_effect, F_sum)
            end
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
    return RunResult(U_steps, U_effect, F_effect)
end

############################################################################################
############################################################################################
"""
    run_threeD_disp(input)

run finite element code for 3D with displacement load
"""
function run_threeD_disp(input::InputStruct)::RunResult
    cv         = input.cell_values
    dh         = input.dh
    ch         = input.ch
    filename   = input.filename
    output_dir = input.output_dir
    
    dof_U = input.dof_U
    dof_F = input.dof_F

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
    F_effect = Vector{Float64}()   
    U_effect = Vector{Float64}() 

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
            apply_zero!(ΔU, ch)
            U .+= ΔU
            #apply!(U, ch)

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

            H1 = "Converged at increment $tot_incr, Δt = $deltaT, total time = $tot_time"
            println(H1);
            if isempty(dof_U) && isempty(dof_F)
                U_effect = []
                F_effect = []
            else
                U_mean = mean(U[dof_U])
                F_sum = sum(F_int[dof_F])
                push!(U_effect, U_mean)
                push!(F_effect, F_sum)
            end
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
    return RunResult( U_steps, U_effect, F_effect )
end
####################

"""
    run_fem(input)

the general fem runner for all cases 2D & 3D 
"""
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