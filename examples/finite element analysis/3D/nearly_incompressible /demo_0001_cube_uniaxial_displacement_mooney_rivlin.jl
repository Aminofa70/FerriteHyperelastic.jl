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

boxDim = [10, 10, 10]
boxEl  = [5, 5, 5]
E, V, F, Fb, Cb = hexbox(boxDim, boxEl)

c1_mod = 10.0
c2_mod = 20.0
ν = 0.4999
μ = 2 * (c1_mod + c2_mod)
K_mod = 2 * μ * (1 + ν) / (3 * (1 - 2 * ν))

displacement_prescribed = 1.0
numSteps = 10

grid = ComodoToFerrite(E, V)

Fb_bottom = Fb[Cb .== 1]
addface!(grid, "bottom", Fb_bottom)

Fb_front = Fb[Cb .== 3]
addface!(grid, "front", Fb_front)

Fb_top = Fb[Cb .== 2]
addface!(grid, "top", Fb_top)

Fb_left = Fb[Cb .== 6]
addface!(grid, "left", Fb_left)


function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "bottom"), (x, t) -> [0.0], [3]))
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "front"),  (x, t) -> [0.0], [2]))
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "left"),   (x, t) -> [0.0], [1]))
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "top"),    (x, t) -> [t],   [3]))
    Ferrite.close!(ch)
    Ferrite.update!(ch, 0.0)
    return ch
end

struct MooneyRivlin
    c1::Float64
    c2::Float64
    K::Float64      # bulk modulus
end


function Ψ(F_arg, Θ, p_field, mp::MooneyRivlin)
    c1, c2, K = mp.c1, mp.c2, mp.K
    J = det(F_arg)
    F̃ = J^(-1/3) * F_arg
    C̃ = tdot(F̃)                   # C̃ = F̃ᵀ F̃

    Ī₁ = tr(C̃)
    Ī₂ = 0.5 * (Ī₁^2 - tr(C̃ ⋅ C̃))

    # --- Deviatoric part ---
    Ψdev = c1 * (Ī₁ - 3.0) + c2 * (Ī₂ - 3.0)

    # --- Volumetric part: U(Θ) = ½ K (ln Θ)² ---
    Ψvol = 0.5 * K * (log(Θ))^2

    # --- Three-field coupling: p(J − Θ)
    Ψcoup = p_field * (J - Θ)

    return Ψdev + Ψvol + Ψcoup
end

function constitutive_driver(F_arg, Θ, p_field, mp::MooneyRivlin)

    # Derivatives w.r.t. F via automatic differentiation
    ∂²Ψ∂F², ∂Ψ∂F = Tensors.hessian(
        y -> Ψ(y, Θ, p_field, mp), F_arg, :all
    )

    # Derivatives w.r.t. Θ — analytical
    #   U(Θ)  = ½ K (ln Θ)²
    #   U'(Θ) = K ln(Θ) / Θ
    #   U''(Θ) = K (1 − ln Θ) / Θ²
    K = mp.K
    lnΘ = log(Θ)
    ∂Ψ∂Θ   = K * lnΘ / Θ - p_field
    ∂²Ψ∂Θ² = K * (1.0 - lnΘ) / (Θ^2)

    return ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂Θ, ∂²Ψ∂Θ²
end


## println(repeat("═", 60))
# ══════════════════════════════════════════════════════════════
# 6. ELEMENT ASSEMBLY
#
#  Three-field: (u, p, Θ)
#
#  Block 1 — u (equilibrium):        
#  Block 2 — p (constraint J = Θ):   
#  Block 3 — Θ (volumetric):         
#
#  
# ══════════════════════════════════════════════════════════════

function assemble_element!(Ke, fe, cell, cvu, cvp, cvΘ,
                           mp, ue, pe, Θe)

    reinit!(cvu, cell)
    reinit!(cvp, cell)
    reinit!(cvΘ, cell)

    fill!(Ke, 0.0)
    fill!(fe, 0.0)

    nu = getnbasefunctions(cvu)
    np = getnbasefunctions(cvp)
    nΘ = getnbasefunctions(cvΘ)

    for qp in 1:getnquadpoints(cvu)
        dΩ = getdetJdV(cvu, qp)

        # --- Current state at quadrature point ---
        ∇u   = function_gradient(cvu, qp, ue)
        F_qp = one(∇u) + ∇u
        p̂    = function_value(cvp, qp, pe)
        Θ̂    = function_value(cvΘ, qp, Θe)

        # --- Constitutive evaluation ---
        ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂Θ, ∂²Ψ∂Θ² =
            constitutive_driver(F_qp, Θ̂, p̂, mp)

        Finv = inv(F_qp)
        J    = det(F_qp)
        # ════════════════════════════════════════════════════
        # Block 1 — u-equation (equilibrium)
        #   Rᵤ[i] = ∫ ∇Nᵢ : P dV,  P = ∂Ψ/∂F
        # ════════════════════════════════════════════════════
        for i in 1:nu
            ∇Nᵢ = shape_gradient(cvu, qp, i)

            # Residual
            fe[BlockIndex((1,),(i,))] += (∇Nᵢ ⊡ ∂Ψ∂F) * dΩ

            # Kuu — material + geometric stiffness
            for j in 1:nu
                ∇Nⱼ = shape_gradient(cvu, qp, j)
                Ke[BlockIndex((1,1),(i,j))] += (∇Nᵢ ⊡ ∂²Ψ∂F² ⊡ ∇Nⱼ) * dΩ
            end

            # δJ_i = J F⁻ᵀ : ∇Nᵢ  (variation of det F)
            δJ_i = J * tr(Finv ⋅ ∇Nᵢ)

            # Kup & Kpu — symmetric coupling
            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)
                Ke[BlockIndex((1,2),(i,j))] += δJ_i * Nⱼᵖ * dΩ
                Ke[BlockIndex((2,1),(j,i))] += δJ_i * Nⱼᵖ * dΩ
            end
        end

        # ════════════════════════════════════════════════════
        # Block 2 — p-equation (constraint J = Θ)
        #   Rₚ[i] = ∫ Nᵢᵖ (J − Θ) dV
        # ════════════════════════════════════════════════════
        for i in 1:np
            Nᵢᵖ = shape_value(cvp, qp, i)

            # Residual
            fe[BlockIndex((2,),(i,))] += Nᵢᵖ * (J - Θ̂) * dΩ

            # KpΘ — off-diagonal coupling
            for j in 1:nΘ
                NⱼΘ = shape_value(cvΘ, qp, j)
                Ke[BlockIndex((2,3),(i,j))] -= Nᵢᵖ * NⱼΘ * dΩ
            end
        end

        # ════════════════════════════════════════════════════
        # Block 3 — Θ-equation (volumetric constitutive)
        #   RΘ[i] = ∫ NᵢΘ (U'(Θ) − p) dV
        #   KΘΘ uses U''(Θ) = K(1 − ln Θ)/Θ²
        # ════════════════════════════════════════════════════
        # 
        for i in 1:nΘ
            NᵢΘ = shape_value(cvΘ, qp, i)

            # Residual
            fe[BlockIndex((3,),(i,))] += NᵢΘ * ∂Ψ∂Θ * dΩ

            # KΘΘ — volumetric tangent
            for j in 1:nΘ
                NⱼΘ = shape_value(cvΘ, qp, j)
                Ke[BlockIndex((3,3),(i,j))] += NᵢΘ * ∂²Ψ∂Θ² * NⱼΘ * dΩ
            end

            # KΘp — (= KpΘᵀ)
            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)
                Ke[BlockIndex((3,2),(i,j))] -= NᵢΘ * Nⱼᵖ * dΩ
            end
        end
    end
end


# ══════════════════════════════════════════════════════════════
# 7. GLOBAL ASSEMBLY
# ══════════════════════════════════════════════════════════════

function assemble_global!(K::SparseMatrixCSC, f,
                          cvu, cvp, cvΘ,
                          dh, mp, w)

    nu = getnbasefunctions(cvu)
    np = getnbasefunctions(cvp)
    nΘ = getnbasefunctions(cvΘ)
    ntot = nu + np + nΘ

    fe = BlockedArray(zeros(ntot), [nu, np, nΘ])
    ke = BlockedArray(zeros(ntot, ntot), [nu, np, nΘ], [nu, np, nΘ])

    assembler = start_assemble(K, f)

    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)

        global_dofs_u = global_dofs[1:nu]
        global_dofs_p = global_dofs[(nu+1):(nu+np)]
        global_dofs_Θ = global_dofs[(nu+np+1):end]

        @assert length(global_dofs) == ntot

        ue = w[global_dofs_u]
        pe = w[global_dofs_p]
        Θe = w[global_dofs_Θ]

        assemble_element!(ke, fe, cell, cvu, cvp, cvΘ,
                          mp, ue, pe, Θe)

        assemble!(assembler, global_dofs, ke, fe)
    end
    return
end


# ══════════════════════════════════════════════════════════════
# 8. FE SETUP HELPERS
# ══════════════════════════════════════════════════════════════

function create_values(ipu, ipp, ipΘ)
    qr  = QuadratureRule{RefHexahedron}(2)
    cvu = CellValues(qr, ipu)
    cvp = CellValues(qr, ipp)
    cvΘ = CellValues(qr, ipΘ)
    return cvu, cvp, cvΘ
end

function create_dofhandler(grid, ipu, ipp, ipΘ)
    dh = DofHandler(grid)
    add!(dh, :u, ipu)
    add!(dh, :p, ipp)
    add!(dh, :Θ, ipΘ)
    close!(dh)
    return dh
end

#Standard Newton iteration for the three-field principle.

function solve(c1, c2, bulk, grid, displacement_prescribed, numSteps;
               NEWTON_TOL = 1e-8, NEWTON_MAXITER = 100)

    # --- Material ---
    mp = MooneyRivlin(c1, c2, bulk)

    # --- Interpolations: Q1/P0/P0 ---
    ipu = Lagrange{RefHexahedron,1}()^3
    ipp = DiscontinuousLagrange{RefHexahedron,0}()
    ipΘ = DiscontinuousLagrange{RefHexahedron,0}()

    # --- FE setup ---
    dh   = create_dofhandler(grid, ipu, ipp, ipΘ)
    dbcs = create_bc(dh)
    cvu, cvp, cvΘ = create_values(ipu, ipp, ipΘ)

    nd = ndofs(dh)

    # --- Storage for visualization ---
    UT         = Vector{Vector{Point{3,Float64}}}(undef, numSteps + 1)
    UT_mag     = Vector{Vector{Float64}}(undef, numSteps + 1)
    ut_mag_max = zeros(Float64, numSteps + 1)

    # --- Solution vectors ---
    un   = zeros(nd)
    u    = zeros(nd)
    Δu   = zeros(nd)
    ΔΔu  = zeros(nd)

    Ksys = allocate_matrix(dh)
    g    = zeros(nd)

    # --- Time parameters ---
    Tf = displacement_prescribed
    Δt = Tf / numSteps

    # --- Initialize Θ DOFs to 1.0 (reference: det F = 1) ---
    Θ_dof_range = dof_range(dh, :Θ)
    for cell in CellIterator(dh)
        gdofs = celldofs(cell)
        for d in gdofs[Θ_dof_range]
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
    disp_pts = [Point{3,Float64}([ux[j], uy[j], uz[j]]) for j in eachindex(ux)]
    UT[1]         = disp_pts
    UT_mag[1]     = norm.(disp_pts)
    ut_mag_max[1] = maximum(UT_mag[1])

    # ══════════════════════════════════════════════════════════
    # TIME STEPPING + NEWTON
    # ══════════════════════════════════════════════════════════
    for step in 1:numSteps
        t = step * Δt
        println("\n=== Time step $step / $numSteps,  t = $t ===")

        Ferrite.update!(dbcs, t)

        fill!(Δu, 0.0)
        newton_itr = 0

        while true
            u .= un .+ Δu
            apply!(u, dbcs)

            assemble_global!(Ksys, g, cvu, cvp, cvΘ, dh, mp, u)

            normg = norm(g[Ferrite.free_dofs(dbcs)])

            if normg < NEWTON_TOL
                println("  Converged in $newton_itr iterations  (‖g‖ = $normg)")
                break
            elseif newton_itr ≥ NEWTON_MAXITER
                error("Newton failed to converge at t = $t")
            end

            apply_zero!(Ksys, g, dbcs)
            ΔΔu .= Ksys \ g
            apply_zero!(ΔΔu, dbcs)

            Δu .-= ΔΔu
            newton_itr += 1
        end

        un .= u

        # --- Post-processing ---
        u_nodes = vec(evaluate_at_grid_nodes(dh, u, :u))
        ux = getindex.(u_nodes, 1)
        uy = getindex.(u_nodes, 2)
        uz = getindex.(u_nodes, 3)
        disp_pts = [Point{3,Float64}([ux[j], uy[j], uz[j]]) for j in eachindex(ux)]

        UT[step + 1]         = disp_pts
        UT_mag[step + 1]     = norm.(disp_pts)
        ut_mag_max[step + 1] = maximum(UT_mag[step + 1])
    end

    return UT, UT_mag, ut_mag_max
end


UT, UT_mag, ut_mag_max = solve(c1_mod, c2_mod, K_mod, grid, displacement_prescribed, numSteps)

numInc = length(UT)

# Create displaced mesh per step
scale = 1.0
VT = [V .+ scale .* UT[i] for i in 1:(numSteps + 1)]
incRange =  0:1:numInc-1

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

hSlider = Slider(fig_disp[2, 1], range=incRange, startvalue= stepStart, linewidth=30)

on(hSlider.value) do stepIndex
    i = stepIndex + 1   # shift to 1-based indexing
    hp[1] = GeometryBasics.Mesh(VT[i], F)
    hp.color = UT_mag[i]
    ax3.title = "Step: $stepIndex"
end

slidercontrol(hSlider, ax3)
display(GLMakie.Screen(), fig_disp)