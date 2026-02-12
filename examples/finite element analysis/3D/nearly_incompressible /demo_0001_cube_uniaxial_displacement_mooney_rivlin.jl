using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics
using IterativeSolvers
using ComodoFerrite
using ComodoFerrite.Ferrite
using SparseArrays
using BlockArrays
## Mesh 
boxDim = [10, 10, 10]
boxEl = [5, 5, 5]
E, V, F, Fb, Cb = hexbox(boxDim, boxEl)
grid = ComodoToFerrite(E, V)

Fb_bottom = Fb[Cb.==1]
addface!(grid , "bottom", Fb_bottom) 

Fb_front = Fb[Cb.==3]  
addface!(grid , "front", Fb_front) 

Fb_top = Fb[Cb.==2] 
addface!(grid , "top", Fb_top)   

Fb_left = Fb[Cb.==6]
addface!(grid , "left", Fb_left)   


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


struct MooneyRivlin
    c1::Float64
    c2::Float64
    K::Float64
end

function Ψ(F, J̄, p, mp::MooneyRivlin)
    c1, c2, K = mp.c1, mp.c2, mp.K

    J  = det(F)
    # Isochoric deformation gradient
    F̄  = J^(-1/3) * F
    # Isochoric right Cauchy-Green: C̄ = F̄ᵀF̄
    C̄  = tdot(F̄)

    # Invariants of C̄
    I1c̄ = tr(C̄)
    I2c̄ = 0.5 * (I1c̄^2 - tr(C̄ ⋅ C̄))

    # Three-field Simo–Taylor–Pister functional:
    #   Π = Ψ_dev(F̄) + Ψ_vol(J̄) + p*(J - J̄)
    Ψdev = c1 * (I1c̄ - 3) + c2 * (I2c̄ - 3)
    Ψvol = 0.5 * K * (log(J̄))^2
    Ψcoup = p * (J - J̄)

    return Ψdev + Ψvol + Ψcoup
end


# ──────────────────────────────────────────────────────────────
# Constitutive driver
# ──────────────────────────────────────────────────────────────
function constitutive_driver(F, J̄, p, mp::MooneyRivlin)

    ∂²Ψ∂F², ∂Ψ∂F = Tensors.hessian(y -> Ψ(y, J̄, p, mp), F, :all)

    ∂²Ψ∂J̄², ∂Ψ∂J̄ = Tensors.hessian(y -> Ψ(F, y, p, mp), J̄, :all)

    return ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂J̄, ∂²Ψ∂J̄²
end

# ──────────────────────────────────────────────────────────────
# Element assembly  (u, p, J̄  three-field Hu–Washizu / STP)
# ──────────────────────────────────────────────────────────────
function assemble_element!(Ke, fe, cell, cvu, cvp, cvJ̄, mp, ue, pe, J̄e)

    reinit!(cvu, cell)
    reinit!(cvp, cell)
    reinit!(cvJ̄, cell)

    fill!(Ke, 0.0)
    fill!(fe, 0.0)

    nu = getnbasefunctions(cvu)
    np = getnbasefunctions(cvp)
    nJ̄ = getnbasefunctions(cvJ̄)

    for qp in 1:getnquadpoints(cvu)
        dΩ = getdetJdV(cvu, qp)

        # Current state at quadrature point
        ∇u = function_gradient(cvu, qp, ue)
        F  = one(∇u) + ∇u
        p̂  = function_value(cvp, qp, pe)
        J̄  = function_value(cvJ̄, qp, J̄e)

        ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂J̄, ∂²Ψ∂J̄² = constitutive_driver(F, J̄, p̂, mp)

        Finv = inv(F)
        J    = det(F)                              

        # ============================================================
        # Block 1 — u-equation  (equilibrium)
        #   Rᵤ[i] = ∫ ∇Nᵢ : P  dV ,   P = ∂Ψ/∂F
        # ============================================================
        for i in 1:nu
            ∇Nᵢ = shape_gradient(cvu, qp, i)

            # --- residual ---
            fe[BlockIndex((1),(i))] += (∇Nᵢ ⊡ ∂Ψ∂F) * dΩ

            # --- Kuu ---
            for j in 1:nu
                ∇Nⱼ = shape_gradient(cvu, qp, j)
                Ke[BlockIndex((1,1),(i,j))] += (∇Nᵢ ⊡ ∂²Ψ∂F² ⊡ ∇Nⱼ) * dΩ
            end

            # δJ_i = J F⁻ᵀ : ∇Nᵢ  (variation of det(F), NOT a J̄ shape fn)
            δJ_i = J * tr(Finv ⋅ ∇Nᵢ)

            # --- Kup  &  Kpu  (symmetric coupling) ---
            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)
                Ke[BlockIndex((1,2),(i,j))] += δJ_i * Nⱼᵖ * dΩ
                Ke[BlockIndex((2,1),(j,i))] += δJ_i * Nⱼᵖ * dΩ
            end
        end

        # ============================================================
        # Block 2 — p-equation  (incompressibility constraint)
        #   Rₚ[i] = ∫ Nᵢᵖ (J − J̄)  dV
        # ============================================================
        for i in 1:np
            Nᵢᵖ = shape_value(cvp, qp, i)

            fe[BlockIndex((2),(i))] += Nᵢᵖ * (J - J̄) * dΩ

            # --- KpJ̄ ---
            for j in 1:nJ̄
                NⱼJ̄ = shape_value(cvJ̄, qp, j)
                Ke[BlockIndex((2,3),(i,j))] -= Nᵢᵖ * NⱼJ̄ * dΩ
            end
        end

        # ============================================================
        # Block 3 — J̄-equation  (volumetric constitutive)
        #   RJ̄[i] = ∫ NᵢJ̄ (∂Ψ_vol/∂J̄ − p)  dV
        # ============================================================
        for i in 1:nJ̄
            NᵢJ̄ = shape_value(cvJ̄, qp, i)

            fe[BlockIndex((3),(i))] += NᵢJ̄ * (∂Ψ∂J̄ - p̂) * dΩ

            # --- KJ̄J̄ ---
            for j in 1:nJ̄
                NⱼJ̄ = shape_value(cvJ̄, qp, j)
                Ke[BlockIndex((3,3),(i,j))] +=
                    NᵢJ̄ * ∂²Ψ∂J̄² * NⱼJ̄ * dΩ
            end

            # --- KJ̄p  (= KpJ̄ᵀ) ---
            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)
                Ke[BlockIndex((3,2),(i,j))] -= NᵢJ̄ * Nⱼᵖ * dΩ
            end
        end
    end
end

function assemble_global!(
        K::SparseMatrixCSC, f,
        cvu::CellValues, cvp::CellValues, cvJ̄::CellValues,
        dh::DofHandler, mp::MooneyRivlin, w
    )

    nu = getnbasefunctions(cvu)
    np = getnbasefunctions(cvp)
    nJ̄ = getnbasefunctions(cvJ̄)
    ntot = nu + np + nJ̄

    # Local blocked arrays: 3 blocks [u, p, J̄]
    fe = BlockedArray(zeros(ntot), [nu, np, nJ̄])
    ke = BlockedArray(zeros(ntot, ntot), [nu, np, nJ̄], [nu, np, nJ̄])

    assembler = start_assemble(K, f)

    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)

        # Extract DOF ranges for each field
        global_dofs_u = global_dofs[1:nu]
        global_dofs_p = global_dofs[(nu + 1):(nu + np)]
        global_dofs_J̄ = global_dofs[(nu + np + 1):end]

        @assert length(global_dofs) == ntot  # sanity check

        # Extract element DOF values from global solution
        ue  = w[global_dofs_u]
        pe  = w[global_dofs_p]
        J̄e  = w[global_dofs_J̄]

        assemble_element!(ke, fe, cell, cvu, cvp, cvJ̄, mp, ue, pe, J̄e)
        assemble!(assembler, global_dofs, ke, fe)
    end
    return
end
 

function create_values(ipu, ipp, ipJ̄)
    qr  = QuadratureRule{RefHexahedron}(2)
    cvu = CellValues(qr, ipu)
    cvp = CellValues(qr, ipp)
    cvJ̄ = CellValues(qr, ipJ̄)
    return cvu, cvp, cvJ̄
end

function create_dofhandler(grid, ipu, ipp, ipJ̄)
    dh = DofHandler(grid)
    add!(dh, :u, ipu)
    add!(dh, :p, ipp)
    add!(dh, :J̄, ipJ̄)             
    close!(dh)
    return dh
end

function solve(c1, c2, bulk, grid, displacement_prescribed, numSteps)
                    

    # --- Material ---
    mp = MooneyRivlin(c1, c2, bulk)                

    # --- Interpolations ---
    # ipu = Lagrange{RefHexahedron, 2}()^3           
    # ipp = Lagrange{RefHexahedron, 1}()              
    # ipJ̄ = Lagrange{RefHexahedron, 1}()


    # Q1/P0/P0 — matches FEBio's hex8 three-field formulation
    ipu = Lagrange{RefHexahedron,1}()^3             # trilinear displacement
    ipp = DiscontinuousLagrange{RefHexahedron,0}()  # element-constant pressure
    ipJ̄ = DiscontinuousLagrange{RefHexahedron,0}()  # element-constant J̄


    # --- FE setup M---
    dh = create_dofhandler(grid, ipu, ipp, ipJ̄)
    dbcs = create_bc(dh)
    cvu, cvp, cvJ̄ = create_values(ipu, ipp, ipJ̄)

    nd = ndofs(dh)

    UT = Vector{Vector{Point{3,Float64}}}(undef, numSteps + 1)
    UT_mag = Vector{Vector{Float64}}(undef, numSteps + 1)
    ut_mag_max = zeros(Float64, numSteps + 1)

    # --- Newton vectors ---
    un = zeros(nd)
    u  = zeros(nd)
    Δu = zeros(nd)
    ΔΔu = zeros(nd)

    Ksys = allocate_matrix(dh)                      
    g = zeros(nd)

    # --- Parameters ---
    NEWTON_TOL = 1e-8
    NEWTON_MAXITER = 100

    Tf = displacement_prescribed
    Δt = Tf / numSteps

    # --- Initial condition ---
    # Initialize J̄ DOFs to 1.0 (reference configuration: det F = 1)
    J̄_dofs = dof_range(dh, :J̄)
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        for d in global_dofs[J̄_dofs]
            un[d] = 1.0
        end
    end

    apply!(un, dbcs)
    u .= un

    # --- Store initial state ---
    u_nodes = vec(evaluate_at_grid_nodes(dh, u, :u))
    ux = getindex.(u_nodes, 1)
    uy = getindex.(u_nodes, 2)
    uz = getindex.(u_nodes, 3)
    disp_points = [Point{3,Float64}([ux[j], uy[j], uz[j]]) for j in eachindex(ux)]
    UT[1] = disp_points
    UT_mag[1] = norm.(disp_points)
    ut_mag_max[1] = maximum(UT_mag[1])

    # --- Time stepping ---
    for step in 1:numSteps
        t = step * Δt
        println("\n=== Time step $step, t = $t ===")

        Ferrite.update!(dbcs, t)

        fill!(Δu, 0.0)
        newton_itr = 0

        while true
            u .= un .+ Δu
            apply!(u, dbcs)
            assemble_global!(Ksys, g, cvu, cvp, cvJ̄, dh, mp, u)

            normg = norm(g[Ferrite.free_dofs(dbcs)])
            if normg < NEWTON_TOL
                println("  Converged in $newton_itr iterations")
                break
            elseif newton_itr ≥ NEWTON_MAXITER
                error("Newton failed to converge at time t = $t")
            end

            apply_zero!(Ksys, g, dbcs)
            ΔΔu .= Ksys \ g

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

        UT[step + 1] = disp_points
        UT_mag[step + 1] = norm.(disp_points)
        ut_mag_max[step + 1] = maximum(UT_mag[step + 1])
    end

    return UT, UT_mag, ut_mag_max
end

c1_mod = 10.0
c2_mod = 20.0
ν = 0.4999
μ = 2*(c1_mod + c2_mod)
K_mod = 2*μ*(1 + ν) / (3*(1 - 2*ν))
K_mod/c1_mod

displacement_prescribed = 2.0
numSteps = 10

UT, UT_mag, ut_mag_max  = solve(c1_mod, c2_mod, K_mod, grid, displacement_prescribed, numSteps)

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