using ComodoFerrite
using ComodoFerrite.Comodo
using ComodoFerrite.Ferrite
using ComodoFerrite.Comodo.GeometryBasics
using ForwardDiff, IterativeSolvers, Roots
using IterativeSolvers, Roots
using ComodoFerrite.Comodo.GLMakie
GLMakie.closeall()

Vec = Ferrite.Vec
## Geometry and Mesh
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "support_1", x -> x[1] ≈ 0.0)
    addfacetset!(grid, "traction", x -> x[1] ≈ Lx)
    return grid
end
## Finite Elemenet Values
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral,order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(2)
    return Ferrite.CellValues(qr, ip), Ferrite.FacetValues(qr_face, ip)
end
## Create Degrees of freedom
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral,1}()^2)
    Ferrite.close!(dh)
    return dh
end

function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "support_1"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.close!(ch)
    return ch
end

struct NeoHooke
    μ::Float64
    λ::Float64
end

function Ψ(C, mp::NeoHooke)
    μ = mp.μ
    λ = mp.λ
    Ic = tr(C)
    J = sqrt(det(C))
    return μ / 2 * (Ic - 3) - μ * log(J) + λ / 2 * (log(J))^2
end

function constitutive_driver(C, mp::NeoHooke)
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, mp), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end;

function assemble_element!(ke, ge, cell, cv, fv, mp, ue, ΓN, tn)
    reinit!(cv, cell)
    fill!(ke, 0.0)
    fill!(ge, 0.0)
    ndofs = getnbasefunctions(cv)

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇u2d = function_gradient(cv, qp, ue)
        F2d = [1.0+∇u2d[1, 1] ∇u2d[1, 2];
            ∇u2d[2, 1] 1.0+∇u2d[2, 2]]

        F = Tensor{2,3,Float64}((
            F2d[1, 1], F2d[1, 2], 0.0,
            F2d[2, 1], F2d[2, 2], 0.0,
            0.0, 0.0, 1.0
        ))
        C = tdot(F) # F' ⋅ F
        # Compute stress and tangent
        S, ∂S∂C = constitutive_driver(C, mp)
        P = F ⋅ S
        I = one(S)
        ∂P∂F = otimesu(I, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(F', I)

        # Loop over test functions
        for i in 1:ndofs
            ∇δui2d = shape_gradient(cv, qp, i)
            ∇δui = Tensor{2,3,Float64}((
                ∇δui2d[1, 1], ∇δui2d[1, 2], 0.0,
                ∇δui2d[2, 1], ∇δui2d[2, 2], 0.0,
                0.0, 0.0, 0.0
            ))

            ge[i] += (∇δui ⊡ P) * dΩ

            ∇δui∂P∂F = ∇δui ⊡ ∂P∂F
            for j in 1:ndofs
                ∇δuj2d = shape_gradient(cv, qp, j)
                ∇δuj = Tensor{2,3,Float64}((
                    ∇δuj2d[1, 1], ∇δuj2d[1, 2], 0.0,
                    ∇δuj2d[2, 1], ∇δuj2d[2, 2], 0.0,
                    0.0, 0.0, 0.0
                ))
                ke[i, j] += (∇δui∂P∂F ⊡ ∇δuj) * dΩ
            end
        end
    end
    for facet in 1:nfacets(cell)
        if (cellid(cell), facet) in ΓN
            reinit!(fv, cell, facet)

            for qp in 1:getnquadpoints(fv)

                ∇uu = function_gradient(fv, qp, ue)
                F = Tensor{2,3,Float64}((1.0 + ∇uu[1, 1], ∇uu[1, 2], 0.0,
                    ∇uu[2, 1], 1.0 + ∇uu[2, 2], 0.0,
                    0.0, 0.0, 1.0))

                J = det(F)
                FinvT = inv(F)'

                C = F' ⋅ F
                Cinv = inv(C)

                N₀_2d = getnormal(fv, qp)
                N₀ = Vec{3}((N₀_2d[1], N₀_2d[2], 0.0))

                m = Cinv ⋅ N₀
                α = N₀ ⋅ m
                sqrtα = sqrt(α)
                Φ = J * sqrtα

                tn_vec = Vec{2}((tn[1], tn[2]))
                T0 = Φ * tn_vec

                dΓ0 = getdetJdV(fv, qp)

                for i in 1:ndofs
                    δui = shape_value(fv, qp, i)
                    ge[i] -= (δui ⋅ T0) * dΓ0

                    for j in 1:ndofs
                        ∇δuj2d = shape_gradient(fv, qp, j)
                        δF = Tensor{2,3,Float64}((
                            ∇δuj2d[1, 1], ∇δuj2d[1, 2], 0.0,
                            ∇δuj2d[2, 1], ∇δuj2d[2, 2], 0.0,
                            0.0, 0.0, 0.0))

                        δJ = J * (FinvT ⊡ δF)
                        δC = δF' ⋅ F + F' ⋅ δF
                        δα = -(m ⋅ (δC ⋅ m))
                        δΦ = sqrtα * δJ + (J / (2 * sqrtα)) * δα
                        δT0 = δΦ * tn_vec

                        ke[i, j] -= (δui ⋅ δT0) * dΓ0
                    end
                end
            end
        end
    end
end;


function assemble_global!(K, g, dh, cv, fv, mp, u, ΓN, tn)
    n = ndofs_per_cell(dh)
    ke = zeros(n, n)
    ge = zeros(n)

    # start_assemble resets K and g
    assembler = start_assemble(K, g)

    # Loop over all cells in the grid
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] # element dofs
        assemble_element!(ke, ge, cell, cv, fv, mp, ue, ΓN, tn)
        assemble!(assembler, global_dofs, ke, ge)
    end
    return
end;

function solve(E, ν, grid, traction_prescribed, numSteps)

    # --- Material ---
    μ = E / (2 * (1 + ν))
    λ = (E * ν) / ((1 + ν) * (1 - 2ν))
    mp = NeoHooke(μ, λ)

    # --- FEM setup ---
    dh = create_dofhandler(grid)
    dbcs = create_bc(dh)
    cv, fv = create_values()

    ΓN = getfacetset(grid, "traction")

    nd = ndofs(dh)

    # --- Storage ---
    UT = Vector{Vector{Point{2,Float64}}}(undef, numSteps + 1)
    UT_mag = Vector{Vector{Float64}}(undef, numSteps + 1)
    ut_mag_max = zeros(numSteps + 1)

    # --- Newton vectors ---
    un = zeros(nd)
    u = zeros(nd)
    Δu = zeros(nd)
    ΔΔu = zeros(nd)

    K = allocate_matrix(dh)
    g = zeros(nd)

    # --- Parameters ---
    NEWTON_TOL = 1e-8
    NEWTON_MAXITER = 100

    Tf = 1.0                    # load factor
    Δt = Tf / numSteps

    # --- Initial condition ---
    apply!(un, dbcs)
    u .= un

    # --- Time stepping ---
    for (step, t) in enumerate(0.0:Δt:Tf)
        println("\n=== Time step $step, load factor = $t ===")

        fill!(Δu, 0.0)
        newton_itr = 0

        # scaled traction vector
        τ = (t * traction_prescribed[1],
            t * traction_prescribed[2])

        while true
            u .= un .+ Δu

            assemble_global!(K, g, dh, cv, fv, mp, u, ΓN, τ)

            normg = norm(g[Ferrite.free_dofs(dbcs)])
            if normg < NEWTON_TOL
                println("  Converged in $newton_itr iterations")
                break
            elseif newton_itr ≥ NEWTON_MAXITER
                error("Newton failed at load factor = $t")
            end

            apply_zero!(K, g, dbcs)

            fill!(ΔΔu, 0.0)
            ΔΔu = K \ g
            apply_zero!(ΔΔu, dbcs)

            Δu .-= ΔΔu
            newton_itr += 1
        end

        un .= u

        # --- Postprocessing ---
        u_nodes = vec(evaluate_at_grid_nodes(dh, u, :u))
        disp_points = [Point{2,Float64}((p[1], p[2])) for p in u_nodes]

        UT[step] = disp_points
        UT_mag[step] = norm.(disp_points)
        ut_mag_max[step] = maximum(UT_mag[step])
    end

    return UT, UT_mag, ut_mag_max
end

Lx, Ly = 1, 1
nx, ny = 10, 10
grid = create_grid(Lx, Ly, nx, ny)
F, V = FerriteToComodo(grid)


E = 10.0
ν = 0.3

traction_prescribed = Vec{2}((5.0, 0.0))
numSteps = 10
UT, UT_mag, ut_mag_max = solve(E, ν, grid, traction_prescribed, numSteps)


# Create displaced mesh per step
scale = 1.0

VV = [Point{2,Float64}(e[1], e[2]) for e in V]
VT = [VV .+ scale .* UT[i] for i in 1:numSteps]

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

# === Visualization setup ===
fig_disp = Figure(size=(1000, 600))
stepStart = 1  # Start at undeformed
ax3 = Axis(fig_disp[1, 1], title="Step: $stepStart")

xlims!(ax3, min_p[1], max_p[1])
ylims!(ax3, min_p[2], max_p[2])
hp = poly!(ax3, GeometryBasics.Mesh(VT[stepStart], F),
    strokewidth=2,
    color=UT_mag[stepStart],
    transparency=false,
    colormap=Reverse(:Spectral),
    colorrange=(0, maximum(ut_mag_max)))

Colorbar(fig_disp[1, 2], hp.plots[1], label="Displacement magnitude [mm]")

incRange = 1:numSteps
hSlider = Slider(fig_disp[2, 1], range=incRange, startvalue=stepStart - 1, linewidth=30)

on(hSlider.value) do stepIndex
    hp[1] = GeometryBasics.Mesh(VT[stepIndex], F)
    hp.color = UT_mag[stepIndex]
    ax3.title = "Step: $stepIndex"
end

slidercontrol(hSlider, ax3)
display(GLMakie.Screen(), fig_disp)
