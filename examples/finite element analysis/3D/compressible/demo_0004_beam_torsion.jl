using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics
using IterativeSolvers
using ComodoFerrite
using ComodoFerrite.Ferrite


## GLMakie setting 
GLMakie.closeall()

## Mesh 
function create_grid(Lx, Ly, Lz, nx, ny, nz)
    left = Ferrite.Vec(0.0, 0.0, 0.0)
    right = Ferrite.Vec(Lx, Ly, Lz)
    grid = generate_grid(Hexahedron, (nx, ny, nz), left, right)
    return grid
end

## Finite Elemenet Values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefHexahedron,order}()^dim
    qr = QuadratureRule{RefHexahedron}(2)
    qr_face = FacetQuadratureRule{RefHexahedron}(1)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

## Create Degrees of freedom
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefHexahedron,1}()^3)
    Ferrite.close!(dh)
    return dh
end

# function rotation(X, t)
#     θ = pi / 3 # 60°
#     x, y, z = X
#     L = 10.0

#     return t * Ferrite.Vec{3}(
#         (
#             L/2 - x + (x - L/2)*cos(θ) + (z - L/2)*sin(θ),
#             0.0,
#             L/2 - z - (x - L/2)*sin(θ) + (z - L/2)*cos(θ),
#         )
#     )
# end

function rotation(X, t; Lx= 10.0, Lz=10.0)
    θ = pi/3 # 60°
    x, y, z = X

    return t * Ferrite.Vec{3}(
        (
            Lx/2 - x + (x - Lx/2)*cos(θ) + (z - Lz/2)*sin(θ),
            0.0,
            Lz/2 - z - (x - Lx/2)*sin(θ) + (z - Lz/2)*cos(θ),
        )
    )
end


function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(dh.grid, "back"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3]) # 
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(dh.grid, "front"), (x, t) -> rotation(x, t), [1, 2, 3]) # 
    add!(ch, dbc)
    Ferrite.close!(ch)
    t = 0.0
    Ferrite.update!(ch, t)
    return ch
end


struct NeoHooke
    μ::Float64
    λ::Float64
end

function Ψ(F, mp::NeoHooke)
    μ = mp.μ
    λ = mp.λ
    Ic = tr(tdot(F))
    J = det(F)
    return μ / 2 * (Ic - 3) - μ * log(J) + λ / 2 * (log(J))^2
end

function constitutive_driver(F, mp::NeoHooke)
    ∂²Ψ∂F², ∂Ψ∂F = Tensors.hessian(y -> Ψ(y,mp), F, :all)
    return ∂Ψ∂F, ∂²Ψ∂F²
end;



function assemble_element!(ke, ge, cell, cv, mp, ue)
    reinit!(cv, cell)
    fill!(ke, 0.0)
    fill!(ge, 0.0)
    ndofs = getnbasefunctions(cv)

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇u = function_gradient(cv, qp, ue)
        F = one(∇u) + ∇u

        ∂Ψ∂F, ∂²Ψ∂F² = constitutive_driver(F, mp)
        P = ∂Ψ∂F

        # Loop over test functions
        for i in 1:ndofs
            ∇δui = shape_gradient(cv, qp, i)

            ge[i] += (∇δui ⊡ P) * dΩ

            ∇δui∂P∂F = ∇δui ⊡ ∂²Ψ∂F²
            for j in 1:ndofs
                ∇δuj = shape_gradient(cv, qp, j)

                ke[i, j] += (∇δui∂P∂F ⊡ ∇δuj) * dΩ
            end
        end
    end
end;

function assemble_global!(K, g, dh, cv, mp, u)
    n = ndofs_per_cell(dh)
    ke = zeros(n, n)
    ge = zeros(n)

    # start_assemble resets K and g
    assembler = start_assemble(K, g)

    # Loop over all cells in the grid
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] # element dofs
        assemble_element!(ke, ge, cell, cv, mp, ue)
        assemble!(assembler, global_dofs, ke, ge)
    end
    return
end;

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


    UT = Vector{Vector{Point{3,Float64}}}(undef, numSteps + 1)
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
    NEWTON_TOL = 1e-8
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
            IterativeSolvers.cg!(ΔΔu, K, g; maxiter=1000)
            apply_zero!(ΔΔu, dbcs)

            Δu .-= ΔΔu
            newton_itr += 1
        end

        un .= u

        # --- Postprocessing ---
        u_nodes = vec(evaluate_at_grid_nodes(dh, u, :u))
        ux = getindex.(u_nodes, 1)
        uy = getindex.(u_nodes, 2)
        uz = getindex.(u_nodes, 3)

        disp_points = [Point{3,Float64}([ux[j], uy[j], uz[j]]) for j in eachindex(ux)]


        UT[step] = disp_points
        UT_mag[step] = norm.(disp_points)
        ut_mag_max[step] = maximum(UT_mag[step])
    end

    return UT, UT_mag, ut_mag_max
end

Lx, Ly, Lz = 10.0, 40.0, 10.0
nx, ny, nz = 5, 10, 5
grid = create_grid(Lx, Ly, Lz, nx, ny, nz)

E, V, F, Fb, Cb = FerriteToComodo(grid, Ferrite.Hexahedron)


E = 10.0
ν = 0.3

displacement_prescribed = 1.0
numSteps = 10
UT, UT_mag, ut_mag_max = solve(E, ν, grid, displacement_prescribed, numSteps)


# Create displaced mesh per step
scale = 1.0
VT = [V .+ scale .* UT[i] for i in 1:numSteps]

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

# === Visualization setup ===
fig_disp = Figure(size=(1000, 600))
stepStart = 1 # Start at undeformed
ax3 = AxisGeom(fig_disp[1, 1], title="Step: $stepStart")


xlims!(ax3, min_p[1], max_p[1])
ylims!(ax3, min_p[2], max_p[2])
zlims!(ax3, min_p[3], max_p[3])

hp = meshplot!(ax3, Fb, VT[stepStart];
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
#display(GLMakie.Screen(), fig_disp)

# frames = vcat(incRange, reverse(incRange[2:end-1]))
# record(fig_disp, "deformation_pingpong.mp4", frames; framerate = 5) do stepIndex
#     hSlider.value[] = stepIndex
# end
