"""
        assemble_cell_3D!( ke_n, fe_int::Vector, cell, cv, input::InputStruct, ue)

This function assembles local internal force and tangent stiffness for 3D
"""
function assemble_cell_3D!( ke_n, fe_int::Vector, cell, cv, input::InputStruct, ue)
  reinit!(cv, cell)
  fill!(ke_n, 0.0)
  fill!(fe_int, 0.0)

  ndofs = getnbasefunctions(cv)
  material = input.material

  for qp in 1:getnquadpoints(cv)
      dΩ = getdetJdV(cv, qp)
      ∇u = function_gradient(cv, qp, ue)
      F = one(∇u) + ∇u

      C = tdot(F)
      S, ∂S∂C = material(C)

      for i in 1:ndofs
          ∇δui = shape_gradient(cv, qp, i)
          δEi = 0.5 * (∇δui' ⋅ F + F' ⋅ ∇δui)

          fe_int[i] += (S ⊡ δEi) * dΩ

          for j in 1:ndofs
              ∇δuj = shape_gradient(cv, qp, j)
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
      assemble_global_3D!(K_nonlinear, F_int, dh, cv, input, u)

This function assembles global internal force and tangent stiffness for 3D
"""
function assemble_global_3D!(K_nonlinear, F_int, dh, cv, input::InputStruct, u)
    n = ndofs_per_cell(dh)
    ke_n = zeros(n, n)
    fe_int = zeros(n)

    # start_assemble resets K and g
    assembler = start_assemble(K_nonlinear, F_int)

    # Loop over all cells in the grid
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] # element dofs
        assemble_cell_3D!(ke_n, fe_int, cell, cv, input, ue)
        assemble!(assembler, global_dofs, ke_n, fe_int)
    end
    return
end;
############################################################################################
############################################################################################