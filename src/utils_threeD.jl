"""
    assemble_cell_3D!(ke_n, fe_int::Vector, cell, cv, input::InputStruct, ue)

Assembles the element stiffness matrix `ke_n` and the internal force vector `fe_int`
for a given finite element `cell` in full 3D.

Arguments
---------
- `ke_n` : Element stiffness matrix to be assembled (modified in place).
- `fe_int::Vector` : Internal force vector for the element (modified in place).
- `cell` : Finite element cell containing nodal connectivity and geometry.
- `cv` : Cell values (e.g., quadrature points, shape function evaluations).
- `input::InputStruct` : Material and model parameters (e.g., constitutive law, elastic modulus).
- `ue` : Vector of nodal displacements for the element.

Notes
-----
- No plane stress/strain assumptions are made: all three displacement components
  and corresponding strains/stresses are active.
- At each quadrature point, the full 3D deformation gradient, stress tensor,
  and material tangent are evaluated.
- Local stiffness and internal force contributions are computed and then assembled
  into the element matrices.
"""
function assemble_cell_3D!(ke_n, fe_int::Vector, cell, cv, input::InputStruct, ue)
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
        I = one(S)
        ∂P∂F = otimesu(I, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(F', I)
        P = F ⋅ S

        for i in 1:ndofs
            ∇δui = shape_gradient(cv, qp, i)

            fe_int[i] += (∇δui ⊡ P) * dΩ

            ∇δui_∂P∂F = ∇δui ⊡ ∂P∂F
            
            for j in 1:ndofs
                ∇δuj = shape_gradient(cv, qp, j)  
                ke_n[i, j] += (∇δui_∂P∂F ⊡ ∇δuj) * dΩ
            end
        end
    end
    return
end 
############################################################################################
############################################################################################
"""
    assemble_global_3D!(K_nonlinear, F_int, dh, cv, input::InputStruct, u)

Assembles the global nonlinear stiffness matrix `K_nonlinear` and the global
internal force vector `F_int` for the entire mesh in full 3D.

Arguments
---------
- `K_nonlinear` : Global stiffness matrix (modified in place).
- `F_int` : Global internal force vector (modified in place).
- `dh` : Discretization handler containing mesh topology, connectivity, and degrees of freedom.
- `cv` : Cell values (e.g., quadrature rules, shape function evaluations for all cells).
- `input::InputStruct` : Material and model parameters (e.g., constitutive law, elastic modulus).
- `u` : Global displacement vector at all degrees of freedom.

Notes
-----
- This function loops over all cells in the 3D mesh and calls the element-level
  routine `assemble_cell_3D!` to compute local stiffness and force contributions.
- No plane stress/strain assumptions are made: all three displacement components
  and corresponding stresses/strains are active.
- The assembled system is consistent with full 3D finite element formulations.
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