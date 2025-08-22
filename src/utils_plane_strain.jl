"""
    assemble_cell_plane_strain!(ke_n, fe_int::Vector, cell, cv, input::InputStruct, ue)

Assembles the element stiffness matrix `ke_n` and the internal force vector `fe_int`
for a given finite element `cell` under plane strain conditions.

Arguments
---------
- `ke_n` : The element stiffness matrix to be assembled (modified in place).
- `fe_int::Vector` : Internal force vector for the element (modified in place).
- `cell` : The finite element cell containing nodal connectivity and geometry.
- `cv` : Cell values (e.g., quadrature, shape function evaluations).
- `input::InputStruct` : Material and model parameters (e.g., constitutive law, elastic modulus).
- `ue` : Vector of nodal displacements for the element.

Notes
-----
This function applies the plane strain assumption, i.e. 
ε_z = 0 but σ_z ≠ 0, typically used for thick structures where 
out-of-plane strains are constrained.

"""
function assemble_cell_plane_strain!(ke_n, fe_int::Vector, cell, cv, input::InputStruct, ue)
    reinit!(cv, cell)
    fill!(ke_n, 0.0)
    fill!(fe_int, 0.0)

    ndofs = getnbasefunctions(cv)    
    material = input.material
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇u2d = function_gradient(cv, qp, ue)
        F = Tensor{2,3,Float64}((
            1.0 + ∇u2d[1, 1], ∇u2d[1, 2], 0.0,
            ∇u2d[2, 1], 1.0 + ∇u2d[2, 2], 0.0,
            0.0, 0.0, 1.0
        ))
        C = tdot(F)
        S, ∂S∂C = material(C)
        I = one(S)
        ∂P∂F = otimesu(I, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(F', I)
        P = F ⋅ S

        for i in 1:ndofs
            ∇δui2d = shape_gradient(cv, qp, i)
            ∇δui = Tensor{2,3,Float64}((
                ∇δui2d[1, 1], ∇δui2d[1, 2], 0.0,
                ∇δui2d[2, 1], ∇δui2d[2, 2], 0.0,
                0.0, 0.0, 0.0
            ))

            fe_int[i] += (∇δui ⊡ P) * dΩ

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

############################################################################################
############################################################################################
"""
    assemble_global_plane_strain!(K_nonlinear, F_int, dh, cv, input::InputStruct, u)

Assembles the global nonlinear stiffness matrix `K_nonlinear` and the global internal
force vector `F_int` for the entire mesh under plane strain conditions.

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
This function loops over all cells in the mesh and calls the element-level assembly
(`assemble_cell_plane_strain!`) to compute local contributions, which are then
assembled into the global stiffness matrix and internal force vector.

The plane strain assumption is applied: ε_z = 0 but σ_z ≠ 0,
suitable for problems where out-of-plane strain is constrained (e.g., thick structures).
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