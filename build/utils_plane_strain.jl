"""
       assemble_cell_plane_strain!( ke_n, fe_int, cell, cv, input, ue)

This function assembles local internal force and tangent stiffness for plane strain
"""
function assemble_cell_plane_strain!( ke_n, fe_int::Vector, cell, cv, input::InputStruct, ue)
    reinit!(cv, cell)
    fill!(ke_n, 0.0)
    fill!(fe_int, 0.0)
    ndofs = getnbasefunctions(cv)
    material = input.material
    
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇uu = function_gradient(cv, qp, ue)

        F = Tensor{2,3,Float64}(( 1.0 + ∇uu[1, 1], ∇uu[1, 2], 0.0, 
                                  ∇uu[2, 1], 1.0 + ∇uu[2, 2], 0.0, 
                                  0.0,        0.0,        1.0 ))
        
        if det(F) <= 0
            error("Jacobian determinant non-positive at qp = $qp")
        end

        C = tdot(F)
        S, ∂S∂C = material(C)

        for i in 1:ndofs
            ∇δui2d = shape_gradient(cv, qp, i)
            ∇δui = Tensor{2,3,Float64}(( ∇δui2d[1,1], ∇δui2d[1,2], 0.0,
                                         ∇δui2d[2,1], ∇δui2d[2,2], 0.0,
                                         0.0,         0.0,        0.0 ))

            δEi = 0.5 * (∇δui' ⋅ F + F' ⋅ ∇δui)
            fe_int[i] += (S ⊡ δEi) * dΩ

            for j in 1:ndofs
                ∇δuj2d = shape_gradient(cv, qp, j)
                ∇δuj = Tensor{2,3,Float64}(( ∇δuj2d[1,1], ∇δuj2d[1,2], 0.0,
                                             ∇δuj2d[2,1], ∇δuj2d[2,2], 0.0,
                                             0.0,         0.0,        0.0 ))

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
       assemble_global_plane_strain!(K_nonlinear, F_int, dh, cv, input, u)

This function assembles global internal force and tangent stiffness for plane strain
"""
function assemble_global_plane_strain!(K_nonlinear, F_int, dh, cv, input::InputStruct, u)
    n = ndofs_per_cell(dh)
    ke_n = zeros(n, n)
    fe_int = zeros(n)

    # start_assemble resets K and g
    assembler = start_assemble(K_nonlinear, F_int)

    # Loop over all cells in the grid
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] # element dofs
        assemble_cell_plane_strain!(ke_n, fe_int, cell, cv, input, ue)
        assemble!(assembler, global_dofs, ke_n, fe_int)
    end
    return
end;
############################################################################################
############################################################################################