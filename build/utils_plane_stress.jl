"""
   solve_lambda3(F2d, input; tol=1e-10, maxit=25)

This functions finds F33 in plane stress 
"""
function solve_lambda3(F2d, input::InputStruct; tol=1e-10, maxit=25)
    
    material = input.material
    # Residual: out-of-plane stress must be zero
    f(λ3) = begin
        F = Tensor{2,3,Float64}((
            F2d[1,1], F2d[1,2], 0.0,
            F2d[2,1], F2d[2,2], 0.0,
            0.0,      0.0,      λ3
        ))
        C = tdot(F)
        S, _ = material(C)
        S[3,3]  # we want S33 == 0
    end

    # Numerical derivative for Newton
    fp(λ3) = (f(λ3 + 1e-8) - f(λ3)) / 1e-8

    # λ3₀ = 1.0
    # Calculate J from 2D deformation gradient using det
    J_2D = det(F2d)
    # Use 1/J_2D as initial guess (for incompressible materials)
    λ3₀ = 1.0 / J_2D

    return find_zero((f, fp), λ3₀, Roots.Newton(); xatol=tol, maxiters=maxit)
end

############################################################################################
############################################################################################
"""
       assemble_cell_plane_stress!(ke_n, fe_int, cell, cv, input, ue)

This function assembles local internal force and tangent stiffness for plane stress
"""
function assemble_cell_plane_stress!(ke_n, fe_int::Vector, cell, cv, input::InputStruct, ue)
    reinit!(cv, cell)
    fill!(ke_n, 0.0)
    fill!(fe_int, 0.0)

    ndofs = getnbasefunctions(cv)
    material = input.material

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇u2d = function_gradient(cv, qp, ue)

        F2d = [1.0 + ∇u2d[1,1]  ∇u2d[1,2];
               ∇u2d[2,1]        1.0 + ∇u2d[2,2]]

        λ3 = solve_lambda3(F2d, input)

        F = Tensor{2,3,Float64}((
            F2d[1,1], F2d[1,2], 0.0,
            F2d[2,1], F2d[2,2], 0.0,
            0.0,      0.0,      λ3
        ))

        C = tdot(F)
        S, ∂S∂C = material(C)

        for i in 1:ndofs
            ∇δui2d = shape_gradient(cv, qp, i)
            ∇δui = Tensor{2,3,Float64}((
                ∇δui2d[1,1], ∇δui2d[1,2], 0.0,
                ∇δui2d[2,1], ∇δui2d[2,2], 0.0,
                0.0,         0.0,        0.0
            ))

            δEi = 0.5 * (∇δui' ⋅ F + F' ⋅ ∇δui)
            fe_int[i] += (S ⊡ δEi) * dΩ

            for j in 1:ndofs
                ∇δuj2d = shape_gradient(cv, qp, j)
                ∇δuj = Tensor{2,3,Float64}((
                    ∇δuj2d[1,1], ∇δuj2d[1,2], 0.0,
                    ∇δuj2d[2,1], ∇δuj2d[2,2], 0.0,
                    0.0,         0.0,        0.0
                ))

                δEj = 0.5 * (∇δuj' ⋅ F + F' ⋅ ∇δuj)

                material_part  = (δEi ⊡ (4 * ∂S∂C) ⊡ δEj)
                geometric_part = (S ⊡ (∇δui' ⋅ ∇δuj))

                ke_n[i, j] += (material_part + geometric_part) * dΩ
            end
        end
    end
    return
end
############################################################################################
############################################################################################
"""
       assemble_global_plane_stress!(K_nonlinear, F_int, dh, cv, input, u)

This function assembles global internal force and tangent stiffness for plane stress
"""
function assemble_global_plane_stress!(K_nonlinear, F_int, dh, cv, input::InputStruct, u)
    n = ndofs_per_cell(dh)
    ke_n = zeros(n, n)
    fe_int = zeros(n)

    # start_assemble resets K and g
    assembler = start_assemble(K_nonlinear, F_int)

    # Loop over all cells in the grid
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] # element dofs
        assemble_cell_plane_stress!(ke_n, fe_int, cell, cv, input, ue)
        assemble!(assembler, global_dofs, ke_n, fe_int)
    end
    return
end;
############################################################################################
############################################################################################