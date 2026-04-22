using FESolvers
using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics
using IterativeSolvers
using ComodoFerrite
using ComodoFerrite.Ferrite
using SparseArrays, LinearAlgebra
using BlockArrays

## GLMakie setting
GLMakie.closeall()


boxDim = [10, 10, 10]
boxEl = [5, 5, 5]
E, V, F, Fb, Cb = hexbox(boxDim, boxEl)
grid = ComodoToFerrite(E, V)

max_displacement = 1.0
c1_mod = 10.0
c2_mod = 20.0
ν = 0.49
μ = 2 * (c1_mod + c2_mod)
K_mod = 2 * μ * (1 + ν) / (3 * (1 - 2 * ν))


Fb_bottom = Fb[Cb .== 1]
addface!(grid, "bottom", Fb_bottom)

Fb_front = Fb[Cb .== 3]
addface!(grid, "front", Fb_front)

Fb_top = Fb[Cb .== 2]
addface!(grid, "top", Fb_top)

Fb_left = Fb[Cb .== 6]
addface!(grid, "left", Fb_left)


struct MooneyRivlin
    c1::Float64
    c2::Float64
    K::Float64
end

function create_values()
    ipu = Lagrange{RefHexahedron, 1}()^3
    ipp = DiscontinuousLagrange{RefHexahedron, 0}()
    ipΘ = DiscontinuousLagrange{RefHexahedron, 0}()
    qr = QuadratureRule{RefHexahedron}(2)
    cvu = CellValues(qr, ipu)
    cvp = CellValues(qr, ipp)
    cvΘ = CellValues(qr, ipΘ)
    return cvu, cvp, cvΘ
end

# DOF HANDLER — three fields (u, p, Θ)
function create_dofhandler(grid)
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefHexahedron, 1}()^3)
    add!(dh, :p, DiscontinuousLagrange{RefHexahedron, 0}())
    add!(dh, :Θ, DiscontinuousLagrange{RefHexahedron, 0}())
    close!(dh)
    return dh
end

function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "bottom"), (x, t) -> [0.0], [3]))
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "front"), (x, t) -> [0.0], [2]))
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "left"), (x, t) -> [0.0], [1]))
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "top"), (x, t) -> [max_displacement * t], [3]))
    Ferrite.close!(ch)
    Ferrite.update!(ch, 0.0)
    return ch
end

#W = U(Θ) + W̃(C̃) + p(J − Θ)
function Ψ(F_arg, Θ, p_field, mp::MooneyRivlin)
    c1, c2, K = mp.c1, mp.c2, mp.K
    J = det(F_arg)

    # Isochoric deformation gradient (eq. 2.28)
    F̃ = J^(-1 / 3) * F_arg
    C̃ = tdot(F̃)

    # Isochoric invariants
    Ī₁ = tr(C̃)
    Ī₂ = 0.5 * (Ī₁^2 - tr(C̃ ⋅ C̃))

    # Deviatoric
    Ψdev = c1 * (Ī₁ - 3.0) + c2 * (Ī₂ - 3.0)
    # Volumetric: U(Θ) = ½ K (ln Θ)²
    Ψvol = 0.5 * K * (log(Θ))^2
    # Coupling: p(J − Θ)
    Ψcoup = p_field * (J - Θ)

    return Ψdev + Ψvol + Ψcoup
end

function constitutive_driver(F_arg, Θ, p_field, mp::MooneyRivlin)
    ∂²Ψ∂F², ∂Ψ∂F = Tensors.hessian(
        y -> Ψ(y, Θ, p_field, mp), F_arg, :all
    )

    K = mp.K
    lnΘ = log(Θ)
    ∂Ψ∂Θ = K * lnΘ / Θ - p_field
    ∂²Ψ∂Θ² = K * (1.0 - lnΘ) / (Θ^2)

    return ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂Θ, ∂²Ψ∂Θ²
end

function assemble_element!(
        Ke, fe, cell, cvu, cvp, cvΘ,
        mp, ue, pe, Θe
    )
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

        ∇u = function_gradient(cvu, qp, ue)
        F_qp = one(∇u) + ∇u
        p̂ = function_value(cvp, qp, pe)
        Θ̂ = function_value(cvΘ, qp, Θe)

        ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂Θ, ∂²Ψ∂Θ² =
            constitutive_driver(F_qp, Θ̂, p̂, mp)

        Finv = inv(F_qp)
        J = det(F_qp)

        # Block 1 — u-equation
        for i in 1:nu
            ∇Nᵢ = shape_gradient(cvu, qp, i)
            fe[BlockIndex((1,), (i,))] += (∇Nᵢ ⊡ ∂Ψ∂F) * dΩ

            for j in 1:nu
                ∇Nⱼ = shape_gradient(cvu, qp, j)
                Ke[BlockIndex((1, 1), (i, j))] += (∇Nᵢ ⊡ ∂²Ψ∂F² ⊡ ∇Nⱼ) * dΩ
            end

            δJ_i = J * tr(Finv ⋅ ∇Nᵢ)
            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)
                Ke[BlockIndex((1, 2), (i, j))] += δJ_i * Nⱼᵖ * dΩ
                Ke[BlockIndex((2, 1), (j, i))] += δJ_i * Nⱼᵖ * dΩ
            end
        end

        # Block 2 — p-equation
        for i in 1:np
            Nᵢᵖ = shape_value(cvp, qp, i)
            fe[BlockIndex((2,), (i,))] += Nᵢᵖ * (J - Θ̂) * dΩ

            for j in 1:nΘ
                NⱼΘ = shape_value(cvΘ, qp, j)
                Ke[BlockIndex((2, 3), (i, j))] -= Nᵢᵖ * NⱼΘ * dΩ
            end
        end

        # Block 3 — Θ-equation
        for i in 1:nΘ
            NᵢΘ = shape_value(cvΘ, qp, i)
            fe[BlockIndex((3,), (i,))] += NᵢΘ * ∂Ψ∂Θ * dΩ

            for j in 1:nΘ
                NⱼΘ = shape_value(cvΘ, qp, j)
                Ke[BlockIndex((3, 3), (i, j))] += NᵢΘ * ∂²Ψ∂Θ² * NⱼΘ * dΩ
            end

            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)
                Ke[BlockIndex((3, 2), (i, j))] -= NᵢΘ * Nⱼᵖ * dΩ
            end
        end
    end
    return
end

function assemble_global!(
        K::SparseMatrixCSC, f,
        cvu, cvp, cvΘ, dh, mp, w
    )
    nu = getnbasefunctions(cvu)
    np = getnbasefunctions(cvp)
    nΘ = getnbasefunctions(cvΘ)
    ntot = nu + np + nΘ

    fe = BlockedArray(zeros(ntot), [nu, np, nΘ])
    ke = BlockedArray(zeros(ntot, ntot), [nu, np, nΘ], [nu, np, nΘ])

    assembler = start_assemble(K, f)

    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        @assert length(global_dofs) == ntot

        ue = w[global_dofs[1:nu]]
        pe = w[global_dofs[(nu + 1):(nu + np)]]
        Θe = w[global_dofs[(nu + np + 1):end]]

        assemble_element!(ke, fe, cell, cvu, cvp, cvΘ, mp, ue, pe, Θe)
        assemble!(assembler, global_dofs, ke, fe)
    end
    return
end

# --- Post-processing storage ---
mutable struct SimoTaylorPost
    solutions::Vector{Vector{Float64}}
    times::Vector{Float64}
end
SimoTaylorPost() = SimoTaylorPost(Vector{Float64}[], Float64[])

# --- Model definition ---
struct SimoTaylorModel{DH, CH, M}
    dh::DH
    ch::CH
    material::M
end

function SimoTaylorModel(mp::MooneyRivlin)
    dh = create_dofhandler(grid)
    ch = create_bc(dh)
    return SimoTaylorModel(dh, ch, mp)
end

# --- Buffer (working arrays) ---
struct SimoTaylorBuffer{CVU, CVP, CVΘ, KT, T}
    cvu::CVU
    cvp::CVP
    cvΘ::CVΘ
    K::KT
    r::Vector{T}
    u::Vector{T}
    time::Vector{T}
end

function build_buffer(model::SimoTaylorModel)
    cvu, cvp, cvΘ = create_values()
    K = allocate_matrix(model.dh)
    nd = ndofs(model.dh)
    r = zeros(nd)
    u = zeros(nd)

    # Initialize Θ DOFs to 1.0
    Θ_dof_range = dof_range(model.dh, :Θ)
    for cell in CellIterator(model.dh)
        gdofs = celldofs(cell)
        for d in gdofs[Θ_dof_range]
            u[d] = 1.0
        end
    end

    return SimoTaylorBuffer(cvu, cvp, cvΘ, K, r, u, [0.0])
end

# --- Problem ---
struct SimoTaylorProblem{PD, PB, PP}
    def::PD
    buf::PB
    post::PP
end

function build_problem(def::SimoTaylorModel)
    return SimoTaylorProblem(def, build_buffer(def), SimoTaylorPost())
end


FESolvers.getunknowns(p::SimoTaylorProblem) = p.buf.u
FESolvers.getresidual(p::SimoTaylorProblem) = p.buf.r
FESolvers.getjacobian(p::SimoTaylorProblem) = p.buf.K

function FESolvers.update_to_next_step!(p::SimoTaylorProblem, time)
    p.buf.time .= time
    Ferrite.update!(p.def.ch, time)
    return apply!(p.buf.u, p.def.ch)
end

function FESolvers.update_problem!(p::SimoTaylorProblem, Δu, _)
    if !isnothing(Δu)
        apply_zero!(Δu, p.def.ch)
        p.buf.u .+= Δu
    end
    assemble_global!(
        p.buf.K, p.buf.r,
        p.buf.cvu, p.buf.cvp, p.buf.cvΘ,
        p.def.dh, p.def.material, p.buf.u
    )
    return apply_zero!(p.buf.K, p.buf.r, p.def.ch)
end

function FESolvers.calculate_convergence_measure(p::SimoTaylorProblem, args...)
    return norm(FESolvers.getresidual(p)[free_dofs(p.def.ch)])
end

function FESolvers.postprocess!(p::SimoTaylorProblem, solver)
    push!(p.post.solutions, copy(p.buf.u))
    push!(p.post.times, p.buf.time[1])
    return println("Time step $(length(p.post.times)) completed, t = $(p.buf.time[1])")
end

FESolvers.handle_converged!(::SimoTaylorProblem) = nothing


mp = MooneyRivlin(c1_mod, c2_mod, K_mod)
def = SimoTaylorModel(mp)
problem = build_problem(def)

solver = QuasiStaticSolver(
    NewtonSolver(; tolerance = 1.0e-6, maxiter = 30),
    AdaptiveTimeStepper(
        0.05, 1.0;
        t_start = 0.0,
        Δt_min = 0.01,
        Δt_max = 0.2
    )
)

solve_problem!(problem, solver)


function solution(prob, grid)
    disp = prob.post.solutions
    numSteps = length(prob.post.times)
    dh = create_dofhandler(grid)

    UT = Vector{Vector{Point{3, Float64}}}(undef, numSteps)
    UT_mag = Vector{Vector{Float64}}(undef, numSteps)
    ut_mag_max = zeros(Float64, numSteps)

    for step in 1:numSteps
        U = disp[step]
        u_nodes = vec(evaluate_at_grid_nodes(dh, U, :u))
        ux = getindex.(u_nodes, 1)
        uy = getindex.(u_nodes, 2)
        uz = getindex.(u_nodes, 3)
        disp_points = [Point{3, Float64}([ux[j], uy[j], uz[j]]) for j in eachindex(ux)]
        UT[step] = disp_points
        UT_mag[step] = norm.(disp_points)
        ut_mag_max[step] = maximum(UT_mag[step])
    end
    return UT, UT_mag, ut_mag_max
end

UT, UT_mag, ut_mag_max = solution(problem, grid)
numSteps = length(problem.post.times)

numInc = length(UT)

# Create displaced mesh per step
scale = 1.0
VT = [V .+ scale .* UT[i] for i in 1:(numSteps)]
incRange = 0:1:(numInc - 1)

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

# === Visualization setup ===
fig_disp = Figure(size = (1000, 600))
stepStart = 1 # Start at undeformed
ax3 = AxisGeom(fig_disp[1, 1], title = "Step: $stepStart")


xlims!(ax3, min_p[1], max_p[1])
ylims!(ax3, min_p[2], max_p[2])
zlims!(ax3, min_p[3], max_p[3])

hp = meshplot!(
    ax3, Fb, VT[stepStart];
    strokewidth = 2,
    color = UT_mag[stepStart],
    transparency = false,
    colormap = Reverse(:Spectral),
    colorrange = (0, maximum(ut_mag_max))
)


Colorbar(fig_disp[1, 2], hp.plots[1], label = "Displacement magnitude [mm]")

hSlider = Slider(fig_disp[2, 1], range = incRange, startvalue = stepStart, linewidth = 30)

on(hSlider.value) do stepIndex
    i = stepIndex + 1   # shift to 1-based indexing
    hp[1] = GeometryBasics.Mesh(VT[i], F)
    hp.color = UT_mag[i]
    ax3.title = "Step: $stepIndex"
end

slidercontrol(hSlider, ax3)
display(GLMakie.Screen(), fig_disp)
