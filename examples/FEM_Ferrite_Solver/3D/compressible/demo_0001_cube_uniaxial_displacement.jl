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
boxDim = [10, 10, 10]
boxEl = [5, 5, 5]
E, V, F, Fb, Cb = hexbox(boxDim, boxEl)
grid = ComodoToFerrite(E, V)

Fb_bottom = Fb[Cb .== 1]
addface!(grid, "bottom", Fb_bottom)

Fb_front = Fb[Cb .== 3]
addface!(grid, "front", Fb_front)

Fb_top = Fb[Cb .== 2]
addface!(grid, "top", Fb_top)

Fb_left = Fb[Cb .== 6]
addface!(grid, "left", Fb_left)


## Finite Elemenet Values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefHexahedron, order}()^dim
    qr = QuadratureRule{RefHexahedron}(2)
    qr_face = FacetQuadratureRule{RefHexahedron}(1)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

## Create Degrees of freedom
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefHexahedron, 1}()^3)
    Ferrite.close!(dh)
    return dh
end

function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(dh.grid, "bottom"), (x, t) -> [0.0], [3]) # bcSupportList_Z
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(dh.grid, "front"), (x, t) -> [0.0], [2]) # bcSupportList_Y
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(dh.grid, "left"), (x, t) -> [0.0], [1]) # bcSupportList_X
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(dh.grid, "top"), (x, t) -> [t], [3]) # bcPrescribeList_Z
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


function assemble_element!(ke, ge, cell, cv, mp, ue)
    reinit!(cv, cell)
    fill!(ke, 0.0)
    fill!(ge, 0.0)
    ndofs = getnbasefunctions(cv)

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇u = function_gradient(cv, qp, ue)
        F = one(∇u) + ∇u

        C = tdot(F) # F' ⋅ F
        # Compute stress and tangent
        S, ∂S∂C = constitutive_driver(C, mp)
        P = F ⋅ S
        I = one(S)
        ∂P∂F = otimesu(I, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(F', I)

        # Loop over test functions
        for i in 1:ndofs
            ∇δui = shape_gradient(cv, qp, i)

            ge[i] += (∇δui ⊡ P) * dΩ

            ∇δui∂P∂F = ∇δui ⊡ ∂P∂F
            for j in 1:ndofs
                ∇δuj = shape_gradient(cv, qp, j)

                ke[i, j] += (∇δui∂P∂F ⊡ ∇δuj) * dΩ
            end
        end
    end
    return
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


function get_top_dofs(dh, grid)
    function dofs_for_component(comp)
        ch_temp = ConstraintHandler(dh)
        add!(ch_temp, Dirichlet(:u, getfacetset(grid, "top"), (x, t) -> [0.0], [comp]))
        close!(ch_temp)
        return ch_temp.prescribed_dofs
    end

    x_dofs = dofs_for_component(1)
    y_dofs = dofs_for_component(2)
    z_dofs = dofs_for_component(3)

    return x_dofs, y_dofs, z_dofs
end

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

    UT = Vector{Vector{Point{3, Float64}}}(undef, numSteps + 1)
    UT_mag = Vector{Vector{Float64}}(undef, numSteps + 1)
    ut_mag_max = zeros(Float64, numSteps + 1)

    Fz = zeros(Float64, numSteps + 1)

    # New curves from solver
    time_curve = zeros(Float64, numSteps + 1)
    disp_curve = zeros(Float64, numSteps + 1)

    # --- Newton vectors ---
    un = zeros(nd)
    u = zeros(nd)
    Δu = zeros(nd)
    ΔΔu = zeros(nd)

    K = allocate_matrix(dh)
    g = zeros(nd)

    NEWTON_TOL = 1.0e-8
    NEWTON_MAXITER = 100

    # --- Initial condition ---
    apply!(un, dbcs)
    u .= un

    x_dofs, y_dofs, z_dofs = get_top_dofs(dh, grid)

    # --- Time stepping ---
    for step in 0:numSteps

        i = step + 1

        # pseudo-time/load factor from 0 to 1
        time = step / numSteps

        # actual prescribed compression displacement
        uz_prescribed = time * displacement_prescribed

        time_curve[i] = time
        disp_curve[i] = uz_prescribed

        println("\n=== Step $step, time = $time, uz = $uz_prescribed ===")

        # Important: update BC with displacement, not physical time
        Ferrite.update!(dbcs, uz_prescribed)

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
                error("Newton failed to converge at step = $step")
            end

            apply_zero!(K, g, dbcs)

            fill!(ΔΔu, 0.0)
            IterativeSolvers.cg!(ΔΔu, K, g; maxiter = 1000)
            apply_zero!(ΔΔu, dbcs)

            Δu .-= ΔΔu
            newton_itr += 1
        end

        un .= u

        # --- Fz reaction force ---
        assemble_global!(K, g, dh, cv, mp, u)
        Fz[i] = sum(g[z_dofs])

        println("  Fz = $(Fz[i])")

        # --- Postprocessing ---
        u_nodes = vec(evaluate_at_grid_nodes(dh, u, :u))

        ux = getindex.(u_nodes, 1)
        uy = getindex.(u_nodes, 2)
        uz = getindex.(u_nodes, 3)

        disp_points = [ Point{3, Float64}([ux[j], uy[j], uz[j]])  for j in eachindex(ux)]

        UT[i] = disp_points
        UT_mag[i] = norm.(disp_points)
        ut_mag_max[i] = maximum(UT_mag[i])
    end

    return UT, UT_mag, ut_mag_max, Fz, time_curve, disp_curve
end


E = 1.0
ν = 0.4

displacement_prescribed = -5.0
numSteps = 10
UT, UT_mag, ut_mag_max, Fz, time_curve, disp_curve = solve(E, ν, grid, displacement_prescribed, numSteps)

Fz_curve = Fz

#########################################################
## Visualization
#########################################################
numInc = length(UT)
scale = 1.0
VT = [V .+ scale .* UT[i] for i in 1:numInc]
min_p = minp([minp(Vi) for Vi in VT])
max_p = maxp([maxp(Vi) for Vi in VT])
incRange = 0:(numInc - 1)
fig = Figure(size = (1600, 800))
stepStart = 1

ax1 = AxisGeom(fig[1, 1], title = "Step: $stepStart", limits = (min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]) )
hp1 = meshplot!(ax1, Fb, VT[stepStart + 1]; strokewidth = 2, color = UT_mag[stepStart + 1], transparency = false, colormap = Reverse(:Spectral), colorrange = (0.0, maximum(ut_mag_max)))
Colorbar(fig[1, 2], hp1.plots[1], label = "Displacement magnitude [mm]")

ax3 = Axis(fig[1, 3], title = "Step: $stepStart", aspect = AxisAspect(1), xlabel = "Time [s]", ylabel = "Force [N]")
lines!(ax3, time_curve, Fz_curve, color = :red, linewidth = 3)

hp3 = scatter!(ax3, Point{2, Float64}(time_curve[stepStart + 1], Fz_curve[stepStart + 1]); markersize = 15, color = :red)
ax4 = Axis(fig[1, 4], title = "Step: $stepStart", aspect = AxisAspect(1), xlabel = "Max displacement [mm]", ylabel = "Force [N]")

lines!(ax4, ut_mag_max, Fz_curve, color = :blue, linewidth = 3)
hp4 = scatter!(ax4, Point{2, Float64}(ut_mag_max[stepStart + 1], Fz_curve[stepStart + 1]); markersize = 15, color = :blue)

hSlider = Slider(fig[2, :], range = incRange, startvalue = stepStart, linewidth = 30)
on(hSlider.value) do stepIndex
    i = stepIndex + 1
    hp1[1] = GeometryBasics.Mesh(VT[i], F)
    hp1.color = UT_mag[i]
    hp3[1] = Point{2, Float64}(time_curve[i], Fz_curve[i])
    hp4[1] = Point{2, Float64}(ut_mag_max[i], Fz_curve[i])
    ax1.title = "Step: $stepIndex"
    ax3.title = "Step: $stepIndex"
    ax4.title = "Step: $stepIndex"
end

slidercontrol(hSlider, ax1)
screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "Neo-Hookean - Reaction Force")
