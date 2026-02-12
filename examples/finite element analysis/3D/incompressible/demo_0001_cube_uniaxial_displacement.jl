using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics
using ComodoFerrite
using ComodoFerrite.Ferrite
using BlockArrays
using Tensors
## GLMakie setting 
GLMakie.closeall()

## Mesh 
function create_grid(Lx, Ly, Lz, nx, ny, nz)
    left = Ferrite.Vec(0.0, 0.0, 0.0)
    right = Ferrite.Vec(Lx, Ly, Lz)
    grid = generate_grid(Ferrite.Tetrahedron, (nx, ny, nz), left, right)
    return grid
end

struct NeoHooke
    μ::Float64
end

function Ψ(F, p, mp::NeoHooke)
    μ = mp.μ
    Ic = tr(tdot(F))
    J  = det(F)
    return μ/2 * (Ic - 3) + p * (J - 1)
end

function constitutive_driver(F, p, mp::NeoHooke)
    ∂²Ψ∂F², ∂Ψ∂F = Tensors.hessian(y -> Ψ(y, p, mp), F, :all)
    ∂²Ψ∂p², ∂Ψ∂p = Tensors.hessian(y -> Ψ(F, y, mp), p, :all)
    ∂²Ψ∂F∂p = Tensors.gradient(q -> Tensors.gradient(y -> Ψ(y, q, mp), F), p)
    return ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂p, ∂²Ψ∂p², ∂²Ψ∂F∂p
end

function assemble_element!(Ke, fe, cell, cellvalues_u, cellvalues_p, mp, ue, pe)
    ublock, pblock = 1, 2
    reinit!(cellvalues_u, cell)
    reinit!(cellvalues_p, cell)
    fill!(Ke, 0.0)
    fill!(fe, 0.0)

    nu = getnbasefunctions(cellvalues_u)
    np = getnbasefunctions(cellvalues_p)

    for qp in 1:getnquadpoints(cellvalues_u)
        dΩ = getdetJdV(cellvalues_u, qp)
        ∇u = function_gradient(cellvalues_u, qp, ue)
        p  = function_value(cellvalues_p, qp, pe)
        F  = one(∇u) + ∇u

        ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂p, _, ∂²Ψ∂F∂p =  constitutive_driver(F, p, mp)

        for i in 1:nu
            ∇δui = shape_gradient(cellvalues_u, qp, i)
            fe[BlockIndex((ublock),(i))] += (∇δui ⊡ ∂Ψ∂F) * dΩ

            for j in 1:nu
                ∇δuj = shape_gradient(cellvalues_u, qp, j)
                Ke[BlockIndex((ublock,ublock),(i,j))] += ((∇δui ⊡ ∂²Ψ∂F²) ⊡ ∇δuj) * dΩ
            end

            for j in 1:np
                δp = shape_value(cellvalues_p, qp, j)
                Ke[BlockIndex((ublock,pblock),(i,j))] += (∂²Ψ∂F∂p ⊡ ∇δui) * δp * dΩ
            end
        end

        for i in 1:np
            δp = shape_value(cellvalues_p, qp, i)
            fe[BlockIndex((pblock),(i))] += δp * ∂Ψ∂p * dΩ

            for j in 1:nu
                ∇δuj = shape_gradient(cellvalues_u, qp, j)
                Ke[BlockIndex((pblock,ublock),(i,j))] += (∇δuj ⊡ ∂²Ψ∂F∂p) * δp * dΩ
            end
        end
    end
end

function assemble_global!(K, f, cellvalues_u, cellvalues_p, dh, mp, w)
    nu = getnbasefunctions(cellvalues_u)
    np = getnbasefunctions(cellvalues_p)

    fe = BlockedArray(zeros(nu+np), [nu,np])
    ke = BlockedArray(zeros(nu+np,nu+np), [nu,np], [nu,np])

    assembler = start_assemble(K, f)
    for cell in CellIterator(dh)
        dofs = celldofs(cell)
        assemble_element!(ke, fe, cell,
                          cellvalues_u, cellvalues_p,
                          mp, w[dofs[1:nu]], w[dofs[nu+1:end]])
        assemble!(assembler, dofs, ke, fe)
    end
end

function create_dofhandler(grid, ipu, ipp)
    dh = DofHandler(grid)
    add!(dh, :u, ipu)
    add!(dh, :p, ipp)
    close!(dh)
    return dh
end

function create_values(interpolation_u, interpolation_p)
    qr = QuadratureRule{RefTetrahedron}(4)
    facet_qr = FacetQuadratureRule{RefTetrahedron}(4)

    cellvalues_u = CellValues(qr, interpolation_u)
    facetvalues_u = FacetValues(facet_qr, interpolation_u)
    cellvalues_p = CellValues(qr, interpolation_p)

    return cellvalues_u, cellvalues_p, facetvalues_u
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

function solve(E, ν, grid, displacement_prescribed, numSteps, interpolation_u, interpolation_p)

    # --- Material ---
    μ = E / (2 * (1 + ν))
    mp = NeoHooke(μ)

    dh = create_dofhandler(grid, interpolation_u, interpolation_p)
    # --- FEM setup ---
    dbcs = create_bc(dh)

    cellvalues_u, cellvalues_p, _ = create_values(interpolation_u, interpolation_p)

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
    f = zeros(nd)

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

            assemble_global!(K, f, cellvalues_u, cellvalues_p, dh, mp, u)

            normg = norm(f[Ferrite.free_dofs(dbcs)])
            if normg < NEWTON_TOL
                println("  Converged in $newton_itr iterations")
                break
            elseif newton_itr ≥ NEWTON_MAXITER
                error("Newton failed to converge at time t = $t")
            end

            apply_zero!(K, f, dbcs)

            fill!(ΔΔu, 0.0)
            ΔΔu = K \ f
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


Lx, Ly, Lz = 10.0, 10.0, 10.0
nx, ny, nz = 5, 5, 5
grid = create_grid(Lx, Ly, Lz, nx, ny, nz)
E, V, F, Fb, Cb = FerriteToComodo(grid)

interpolation_u = Lagrange{RefTetrahedron,2}()^3
interpolation_p = Lagrange{RefTetrahedron,1}()

E = 10.0
ν = 0.5

displacement_prescribed = 5.0
numSteps = 10
UT, UT_mag, ut_mag_max = solve(E, ν, grid, displacement_prescribed, numSteps, interpolation_u, interpolation_p)


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
               strokewidth = 2,
               color = UT_mag[stepStart], 
               transparency = false, 
               colormap = Reverse(:Spectral), 
               colorrange = (0, maximum(ut_mag_max)))


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

