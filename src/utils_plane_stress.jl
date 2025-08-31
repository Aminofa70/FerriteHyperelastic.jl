"""
    solve_lambda3(F2d, input::InputStruct; tol=1e-10, maxit=25)

Solves for the out-of-plane stretch `λ₃` required to satisfy the
plane stress condition (σ₃₃ = 0), given the 2D deformation gradient `F2d`.

Arguments
---------
- `F2d` : 2×2 in-plane deformation gradient matrix.
- `input::InputStruct` : Material and model parameters (e.g., constitutive law, compressibility).
- `tol` : Convergence tolerance for the iterative procedure (default: 1e-10).
- `maxit` : Maximum number of iterations allowed (default: 25).

Returns
-------
- `λ₃` : The out-of-plane stretch ensuring the plane stress condition σ₃₃ = 0.

Notes
-----
This function solves a nonlinear scalar equation for `λ₃`, typically using
a Newton–Raphson iteration. The residual is based on the out-of-plane stress
component, and iterations stop once |σ₃₃| falls below `tol` or after `maxit`
iterations.
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
    assemble_cell_plane_stress!(ke_n, fe_int::Vector, cell, cv, input::InputStruct, ue)

Assembles the element stiffness matrix `ke_n` and the internal force vector `fe_int`
for a given finite element `cell` under plane stress conditions.

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
- The plane stress condition (σ₃₃ = 0) is enforced by solving for the out-of-plane stretch `λ₃`
  at each quadrature point using `solve_lambda3`.  
- Local stiffness and force contributions are computed from the 2D deformation state,
  adjusted with the solved `λ₃`, and then assembled into the element matrices.
"""
function assemble_cell_plane_stress!(
    ke_n,
    fe_int::Vector,
    cell,
    cv,
    input::InputStruct,
    ue
)
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
    assemble_global_plane_stress!(K_nonlinear, F_int, dh, cv, input::InputStruct, u)

Assembles the global nonlinear stiffness matrix `K_nonlinear` and the global
internal force vector `F_int` for the entire mesh under plane stress conditions.

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
- This function loops over all cells in the mesh and calls the element-level assembly
  (`assemble_cell_plane_stress!`) to compute local stiffness and force contributions.
- The plane stress condition (σ₃₃ = 0) is enforced locally within each element
  by solving for the out-of-plane stretch `λ₃` using `solve_lambda3`.
- The assembled global system is consistent with 2D plane stress formulations,
  where out-of-plane stresses vanish but out-of-plane strains may be nonzero.
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