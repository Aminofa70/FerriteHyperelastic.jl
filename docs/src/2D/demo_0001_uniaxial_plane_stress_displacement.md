# Finite Element for uniaxial plane stress subjected to prescribed displacement

### Dependencies

````julia
using ComodoFerrite
using ComodoFerrite.Comodo
using ComodoFerrite.Ferrite
using ComodoFerrite.Comodo.GeometryBasics
using ForwardDiff, IterativeSolvers, Roots
using IterativeSolvers, Roots
using ComodoFerrite.Comodo.GLMakie
````

### Initialization of the plots

````julia
GLMakie.closeall()
````

### Creating a Rectangular Finite Element Grid with Boundary Node Sets

````julia
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly)),
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "support_1", x -> x[1] ≈ 0.0)
    addnodeset!(grid, "load", x -> x[1] ≈ Lx)
    return grid
end
````

###  Finite Element Shape Functions and Quadrature Setup

````julia
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral, order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(2)
    return Ferrite.CellValues(qr, ip), Ferrite.FacetValues(qr_face, ip)
end
````

### Degree of Freedom (DoF) Handler Setup

````julia
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()^2)
    Ferrite.close!(dh)
    return dh
end
````

### Boundary Condition (Dirichlet) Setup

````julia
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "support_1"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "load"), (x, t) -> [t], [1]))
    Ferrite.close!(ch)
    t = 0.0
    Ferrite.update!(ch, t)
    return ch
end
````

### Neo-Hookean Material Model and Strain Energy Function

````julia
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
````

### Constitutive Driver: Stress and Tangent Computation

````julia
function constitutive_driver(C, mp::NeoHooke)
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, mp), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end;
````

### Solving Out-of-Plane Stretch (λ₃) via Newton Method

````julia
function solve_lambda3(F2d, mp; tol = 1.0e-10, maxit = 25)
    J2D = det(F2d)
    λ3₀ = inv(J2D)

    function residual(λ3::T) where {T <: Real}
        Z = zero(T)

        F = Tensor{2, 3, T}(
            (
                T(F2d[1, 1]), T(F2d[1, 2]), Z,
                T(F2d[2, 1]), T(F2d[2, 2]), Z,
                Z, Z, λ3,
            )
        )
        if det(F) <= 0
            error("Jacobian determinant non-positive at qp = $qp")
        end

        C = tdot(F)
        S, _ = constitutive_driver(C, mp)
        return S[3, 3]
    end

    jacobian(λ3) = ForwardDiff.derivative(residual, λ3)

    return find_zero(
        (residual, jacobian),
        λ3₀,
        Roots.Newton();
        xatol = tol,
        maxiters = maxit
    )
end
````

### Element Residual and Tangent Stiffness Assembly

````julia
function assemble_element!(ke, ge, cell, cv, mp, ue)
    reinit!(cv, cell)
    fill!(ke, 0.0)
    fill!(ge, 0.0)
    ndofs = getnbasefunctions(cv)

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇u2d = function_gradient(cv, qp, ue)
        F2d = [
            1.0 + ∇u2d[1, 1] ∇u2d[1, 2];
            ∇u2d[2, 1] 1.0 + ∇u2d[2, 2]
        ]


        λ3 = solve_lambda3(F2d, mp)

        F = Tensor{2, 3, Float64}(
            (
                F2d[1, 1], F2d[1, 2], 0.0,
                F2d[2, 1], F2d[2, 2], 0.0,
                0.0, 0.0, λ3,
            )
        )
        C = tdot(F) ## F' ⋅ F
        # Compute stress and tangent
        S, ∂S∂C = constitutive_driver(C, mp)
        P = F ⋅ S
        I = one(S)
        ∂P∂F = otimesu(I, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(F', I)

        # Loop over test functions
        for i in 1:ndofs
            ∇δui2d = shape_gradient(cv, qp, i)
            ∇δui = Tensor{2, 3, Float64}(
                (
                    ∇δui2d[1, 1], ∇δui2d[1, 2], 0.0,
                    ∇δui2d[2, 1], ∇δui2d[2, 2], 0.0,
                    0.0, 0.0, 0.0,
                )
            )

            ge[i] += (∇δui ⊡ P) * dΩ

            ∇δui∂P∂F = ∇δui ⊡ ∂P∂F
            for j in 1:ndofs
                ∇δuj2d = shape_gradient(cv, qp, j)
                ∇δuj = Tensor{2, 3, Float64}(
                    (
                        ∇δuj2d[1, 1], ∇δuj2d[1, 2], 0.0,
                        ∇δuj2d[2, 1], ∇δuj2d[2, 2], 0.0,
                        0.0, 0.0, 0.0,
                    )
                )
                ke[i, j] += (∇δui∂P∂F ⊡ ∇δuj) * dΩ
            end
        end
    end
    return
end;
````

### Global Residual and Stiffness Matrix Assembly

````julia
function assemble_global!(K, g, dh, cv, mp, u)
    n = ndofs_per_cell(dh)
    ke = zeros(n, n)
    ge = zeros(n)

    # start_assemble resets K and g
    assembler = start_assemble(K, g)

    # Loop over all cells in the grid
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] ## element dofs
        assemble_element!(ke, ge, cell, cv, mp, ue)
        assemble!(assembler, global_dofs, ke, ge)
    end
    return
end;
````

### Nonlinear Finite Element Solver with Load Stepping

````julia
function solve(E, ν, grid, displacement_prescribed, numSteps)

    # --- Material ---
    μ = E / (2 * (1 + ν))
    λ = (E * ν) / ((1 + ν) * (1 - 2ν))
    mp = NeoHooke(μ, λ)

    # --- FEM setup ---
    dh = create_dofhandler(grid)
    dbcs = create_bc(dh)
    cv, _ = create_values()


    nd = ndofs(dh)


    UT = Vector{Vector{Point{2, Float64}}}(undef, numSteps + 1)
    UT_mag = Vector{Vector{Float64}}(undef, numSteps + 1)
    ut_mag_max = zeros(Float64, numSteps + 1)

    # --- Newton vectors ---
    un = zeros(nd)
    u = zeros(nd)
    Δu = zeros(nd)
    ΔΔu = zeros(nd)

    K = allocate_matrix(dh)
    g = zeros(nd)

    # --- Parameters ---
    NEWTON_TOL = 1.0e-8
    NEWTON_MAXITER = 100

    Tf = displacement_prescribed
    Δt = Tf / numSteps

    # --- Initial condition ---
    apply!(un, dbcs)
    u .= un

    # --- Time stepping (UNCHANGED) ---
    for (step, t) in enumerate(0.0:Δt:Tf)
        println("\n=== Time step $step, t = $t ===")

        Ferrite.update!(dbcs, t)

        fill!(Δu, 0.0)
        newton_itr = 0

        while true
            u .= un .+ Δu
            apply!(u, dbcs)

            assemble_global!(K, g, dh, cv, mp, u)

            normg = norm(g[Ferrite.free_dofs(dbcs)])
            if normg < NEWTON_TOL
                println("  Converged in $newton_itr iterations")
                break
            elseif newton_itr ≥ NEWTON_MAXITER
                error("Newton failed to converge at time t = $t")
            end

            apply_zero!(K, g, dbcs)

            fill!(ΔΔu, 0.0)
            IterativeSolvers.cg!(ΔΔu, K, g; maxiter = 1000)
            apply_zero!(ΔΔu, dbcs)

            Δu .-= ΔΔu
            newton_itr += 1
        end

        un .= u

        # --- Postprocessing ---
        u_nodes = vec(evaluate_at_grid_nodes(dh, u, :u))
        ux = getindex.(u_nodes, 1)
        uy = getindex.(u_nodes, 2)

        disp_points = [Point{2, Float64}((ux[i], uy[i])) for i in eachindex(ux)]

        UT[step] = disp_points
        UT_mag[step] = norm.(disp_points)
        ut_mag_max[step] = maximum(UT_mag[step])
    end

    return UT, UT_mag, ut_mag_max
end
````

### Grid Creation and Conversion to Comodo Format

````julia
Lx, Ly = 1, 1
nx, ny = 10, 10
grid = create_grid(Lx, Ly, nx, ny)
F, V = FerriteToComodo(grid)
````

### Material Properties Definition (Young’s Modulus and Poisson’s Ratio)

````julia
E = 10.0
ν = 0.3
````

### Running the Simulation and Computing Displacement Results

````julia
displacement_prescribed = 1.0
numSteps = 10
UT, UT_mag, ut_mag_max = solve(E, ν, grid, displacement_prescribed, numSteps)
````

### Postprocessing Deformed Mesh and Animation Setup

````julia
numInc = length(UT)          ## 11 steps: 0 → 10
scale = 1.0

# Convert V to 2D points
VV = [Point{2, Float64}(e[1], e[2]) for e in V]

# Use all steps (1-based indexing in Julia)
VT = [VV .+ scale .* UT[i] for i in 1:numInc]

# Slider from 0 → numInc-1
incRange = 0:(numInc - 1)

min_p = minp([minp(VT[i]) for i in 1:numInc])
max_p = maxp([maxp(VT[i]) for i in 1:numInc])
````

### Interactive Visualization of Deformation with Slider Control

````julia
fig_disp = Figure(size = (1000, 600))
stepStart = 1
ax3 = Axis(fig_disp[1, 1], title = "Step: $stepStart")

xlims!(ax3, min_p[1], max_p[1])
ylims!(ax3, min_p[2], max_p[2])

# Initial mesh (step 0)
hp = poly!(
    ax3, GeometryBasics.Mesh(VT[stepStart + 1], F),
    strokewidth = 2,
    color = UT_mag[stepStart + 1],
    transparency = false,
    colormap = Reverse(:Spectral),
    colorrange = (0, maximum(ut_mag_max))
)

Colorbar(fig_disp[1, 2], hp.plots[1], label = "Displacement magnitude [mm]")

hSlider = Slider(fig_disp[2, 1], range = incRange, startvalue = stepStart, linewidth = 30)

on(hSlider.value) do stepIndex
    i = stepIndex + 1  ## shift to 1-based array index
    hp[1] = GeometryBasics.Mesh(VT[i], F)
    hp.color = UT_mag[i]
    ax3.title = "Step: $stepIndex"
end

slidercontrol(hSlider, ax3)
display(GLMakie.Screen(), fig_disp)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

