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

const FVec = Ferrite.Vec

## GLMakie setting
GLMakie.closeall()

function create_cook_grid(nx, ny)
    corners = [
        FVec{2}((0.0, 0.0)),
        FVec{2}((48.0, 44.0)),
        FVec{2}((48.0, 60.0)),
        FVec{2}((0.0, 44.0)),
    ]
    grid = generate_grid(Quadrilateral, (nx, ny), corners)
    addfacetset!(grid, "clamped", x -> norm(x[1]) ≈ 0.0)
    addfacetset!(grid, "traction", x -> norm(x[1]) ≈ 48.0)
    addnodeset!(grid, "pointA", x -> x[1] ≈ 48.0 && x[2] ≈ 60.0)
    return grid
end

grid = create_cook_grid(32, 32)


μ  = 80.1938   # MPa
ν  = 0.4999
κ  = 2.0 * μ * (1.0 + ν) / (3.0 * (1.0 - 2.0 * ν))

traction_value = 24.0  # N/mm²

struct NeoHookean
    μ::Float64
    κ::Float64
end


function create_values()
    ipu = Lagrange{RefQuadrilateral, 1}()^2
    ipp = DiscontinuousLagrange{RefQuadrilateral, 0}()
    ipΘ = DiscontinuousLagrange{RefQuadrilateral, 0}()
    qr  = QuadratureRule{RefQuadrilateral}(2)
    qr_face = FacetQuadratureRule{RefQuadrilateral}(2)
    cvu = CellValues(qr, ipu)
    cvp = CellValues(qr, ipp)
    cvΘ = CellValues(qr, ipΘ)
    fv  = FacetValues(qr_face, ipu)
    return cvu, cvp, cvΘ, fv
end


function create_dofhandler(grid)
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral, 1}()^2)
    add!(dh, :p, DiscontinuousLagrange{RefQuadrilateral, 0}())
    add!(dh, :Θ, DiscontinuousLagrange{RefQuadrilateral, 0}())
    close!(dh)
    return dh
end

function create_bc(dh)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "clamped"), (x, t) -> [0.0, 0.0], [1, 2]))
    close!(ch)
    return ch
end

function Ψ(F3d::Tensor{2,3}, Θ::Number, p_field::Number, mp::NeoHookean)
    μ_val, κ_val = mp.μ, mp.κ
    J = det(F3d)

    # Isochoric deformation gradient
    F̃ = J^(-1 / 3) * F3d
    C̃ = tdot(F̃)
    Ī₁ = tr(C̃)

    # ψ_iso = μ/2 (Ī₁ - 3)
    Ψiso = μ_val / 2.0 * (Ī₁ - 3.0)
    # U(Θ) = κ/4 (Θ² - 1 - 2 ln Θ)
    Ψvol = κ_val / 4.0 * (Θ^2 - 1.0 - 2.0 * log(Θ))
    # p(J - Θ)
    Ψcoup = p_field * (J - Θ)

    return Ψiso + Ψvol + Ψcoup
end


function constitutive_driver(F3d, Θ, p_field, mp::NeoHookean)
    ∂²Ψ∂F², ∂Ψ∂F = Tensors.hessian(
        y -> Ψ(y, Θ, p_field, mp), F3d, :all
    )

    κ_val = mp.κ
    ∂Ψ∂Θ   = κ_val / 2.0 * (Θ - 1.0 / Θ) - p_field
    ∂²Ψ∂Θ² = κ_val / 2.0 * (1.0 + 1.0 / (Θ^2))

    return ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂Θ, ∂²Ψ∂Θ²
end

function assemble_element!(Ke, fe, cell, cvu, cvp, cvΘ, fv,
                           mp, ue, pe, Θe, ΓN, tn, time)
    reinit!(cvu, cell)
    reinit!(cvp, cell)
    reinit!(cvΘ, cell)

    fill!(Ke, 0.0)
    fill!(fe, 0.0)

    nu = getnbasefunctions(cvu)
    np = getnbasefunctions(cvp)
    nΘ = getnbasefunctions(cvΘ)

    # ---- Volume integration ----
    for qp in 1:getnquadpoints(cvu)
        dΩ = getdetJdV(cvu, qp)

        # 2D gradient → 3D (plane strain: F33 = 1)
        ∇u2d = function_gradient(cvu, qp, ue)
        F_qp = Tensor{2,3,Float64}((
            1.0 + ∇u2d[1,1], ∇u2d[1,2], 0.0,
            ∇u2d[2,1], 1.0 + ∇u2d[2,2], 0.0,
            0.0, 0.0, 1.0
        ))

        p̂ = function_value(cvp, qp, pe)
        Θ̂ = function_value(cvΘ, qp, Θe)

        ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂Θ, ∂²Ψ∂Θ² =
            constitutive_driver(F_qp, Θ̂, p̂, mp)

        Finv = inv(F_qp)
        J    = det(F_qp)

        # Block 1 — u-equation
        for i in 1:nu
            ∇Nᵢ2d = shape_gradient(cvu, qp, i)
            ∇Nᵢ = Tensor{2,3,Float64}((
                ∇Nᵢ2d[1,1], ∇Nᵢ2d[1,2], 0.0,
                ∇Nᵢ2d[2,1], ∇Nᵢ2d[2,2], 0.0,
                0.0, 0.0, 0.0
            ))

            fe[BlockIndex((1,),(i,))] += (∇Nᵢ ⊡ ∂Ψ∂F) * dΩ

            for j in 1:nu
                ∇Nⱼ2d = shape_gradient(cvu, qp, j)
                ∇Nⱼ = Tensor{2,3,Float64}((
                    ∇Nⱼ2d[1,1], ∇Nⱼ2d[1,2], 0.0,
                    ∇Nⱼ2d[2,1], ∇Nⱼ2d[2,2], 0.0,
                    0.0, 0.0, 0.0
                ))
                Ke[BlockIndex((1,1),(i,j))] += (∇Nᵢ ⊡ ∂²Ψ∂F² ⊡ ∇Nⱼ) * dΩ
            end

            δJ_i = J * tr(Finv ⋅ ∇Nᵢ)
            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)
                Ke[BlockIndex((1,2),(i,j))] += δJ_i * Nⱼᵖ * dΩ
                Ke[BlockIndex((2,1),(j,i))] += δJ_i * Nⱼᵖ * dΩ
            end
        end

        # Block 2 — p-equation
        for i in 1:np
            Nᵢᵖ = shape_value(cvp, qp, i)
            fe[BlockIndex((2,),(i,))] += Nᵢᵖ * (J - Θ̂) * dΩ

            for j in 1:nΘ
                NⱼΘ = shape_value(cvΘ, qp, j)
                Ke[BlockIndex((2,3),(i,j))] -= Nᵢᵖ * NⱼΘ * dΩ
            end
        end

        # Block 3 — Θ-equation
        for i in 1:nΘ
            NᵢΘ = shape_value(cvΘ, qp, i)
            fe[BlockIndex((3,),(i,))] += NᵢΘ * ∂Ψ∂Θ * dΩ

            for j in 1:nΘ
                NⱼΘ = shape_value(cvΘ, qp, j)
                Ke[BlockIndex((3,3),(i,j))] += NᵢΘ * ∂²Ψ∂Θ² * NⱼΘ * dΩ
            end

            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)
                Ke[BlockIndex((3,2),(i,j))] -= NᵢΘ * Nⱼᵖ * dΩ
            end
        end
    end

    # ---- Dead load traction on Neumann boundary ----
    tn_current = time * tn
    for facet in 1:nfacets(cell)
        if FacetIndex(cellid(cell), facet) in ΓN
            reinit!(fv, cell, facet)
            for qp in 1:getnquadpoints(fv)
                dΓ0 = getdetJdV(fv, qp)
                for i in 1:nu
                    δui = shape_value(fv, qp, i)
                    fe[BlockIndex((1,),(i,))] -= (δui ⋅ tn_current) * dΓ0
                end
            end
        end
    end
end


function assemble_global!(K::SparseMatrixCSC, f,
                          cvu, cvp, cvΘ, fv, dh, mp, w, ΓN, tn, time)
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
        pe = w[global_dofs[(nu+1):(nu+np)]]
        Θe = w[global_dofs[(nu+np+1):end]]

        assemble_element!(ke, fe, cell, cvu, cvp, cvΘ, fv,
                          mp, ue, pe, Θe, ΓN, tn, time)
        assemble!(assembler, global_dofs, ke, fe)
    end
    return
end

# ============================================================
# Post-processing storage
# ============================================================
mutable struct CookPost
    solutions::Vector{Vector{Float64}}
    times::Vector{Float64}
end
CookPost() = CookPost(Vector{Float64}[], Float64[])

# ============================================================
# Model definition
# ============================================================
struct CookModel{DH,CH,FS}
    dh::DH
    ch::CH
    material::NeoHookean
    ΓN::FS
    tn::FVec{2,Float64}
end

# ============================================================
# Buffer
# ============================================================
struct CookBuffer{CVU,CVP,CVΘ,FV,KT,T}
    cvu::CVU
    cvp::CVP
    cvΘ::CVΘ
    fv::FV
    K::KT
    r::Vector{T}
    u::Vector{T}
    time::Vector{T}
end

function build_buffer(model::CookModel)
    cvu, cvp, cvΘ, fv = create_values()
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

    return CookBuffer(cvu, cvp, cvΘ, fv, K, r, u, [0.0])
end


struct CookProblem{PD,PB,PP}
    def::PD
    buf::PB
    post::PP
end

function build_problem(def::CookModel)
    CookProblem(def, build_buffer(def), CookPost())
end

FESolvers.getunknowns(p::CookProblem) = p.buf.u
FESolvers.getresidual(p::CookProblem) = p.buf.r
FESolvers.getjacobian(p::CookProblem) = p.buf.K

function FESolvers.update_to_next_step!(p::CookProblem, time)
    p.buf.time .= time
    Ferrite.update!(p.def.ch, time)
    apply!(p.buf.u, p.def.ch)
end

function FESolvers.update_problem!(p::CookProblem, Δu, _)
    if !isnothing(Δu)
        apply_zero!(Δu, p.def.ch)
        p.buf.u .+= Δu
    end
    assemble_global!(p.buf.K, p.buf.r,
                     p.buf.cvu, p.buf.cvp, p.buf.cvΘ, p.buf.fv,
                     p.def.dh, p.def.material, p.buf.u,
                     p.def.ΓN, p.def.tn, p.buf.time[1])
    apply_zero!(p.buf.K, p.buf.r, p.def.ch)
end

function FESolvers.calculate_convergence_measure(p::CookProblem, args...)
    norm(FESolvers.getresidual(p)[free_dofs(p.def.ch)])
end

function FESolvers.postprocess!(p::CookProblem, solver)
    push!(p.post.solutions, copy(p.buf.u))
    push!(p.post.times, p.buf.time[1])
    println("Time step $(length(p.post.times)) completed, t = $(p.buf.time[1])")
end

FESolvers.handle_converged!(::CookProblem) = nothing

mp = NeoHookean(μ, κ)

dh = create_dofhandler(grid)
ch = create_bc(dh)

ΓN = getfacetset(grid, "traction")
tn = FVec{2}((0.0, traction_value))

def = CookModel(dh, ch, mp, ΓN, tn)
problem = build_problem(def)

solver = QuasiStaticSolver(
    NewtonSolver(; tolerance=1e-6, maxiter=30),
    FixedTimeStepper(; num_steps=10, Δt=0.1, t_start=0.0)
)

solve_problem!(problem, solver)

# Extract tip displacement at point A (48, 60)

dh_post = create_dofhandler(grid)

ch_A = ConstraintHandler(dh_post)
add!(ch_A, Dirichlet(:u, getnodeset(grid, "pointA"), (x, t) -> [0.0, 0.0], [1, 2]))
close!(ch_A)
pointA_dofs = ch_A.prescribed_dofs

numSteps = length(problem.post.times)
u1_A = zeros(numSteps)
u2_A = zeros(numSteps)

for step in 1:numSteps
    u = problem.post.solutions[step]
    u1_A[step] = u[pointA_dofs[1]]
    u2_A[step] = u[pointA_dofs[2]]
end

println("\n=== Tip Displacement at Point A (48, 60) ===")
for step in 1:numSteps
    println("Step $step, t = $(problem.post.times[step]): u₁ = $(u1_A[step]),  u₂ = $(u2_A[step])")
end
println("\nFinal tip vertical displacement u₂ = ", u2_A[end])

F, V   = FerriteToComodo(grid)

function extract_solution(prob, grid)
    disp = prob.post.solutions
    numSteps = length(prob.post.times)
    dh = create_dofhandler(grid)

    UT     = Vector{Vector{Point{2,Float64}}}(undef, numSteps)
    UT_mag = Vector{Vector{Float64}}(undef, numSteps)
    ut_mag_max = zeros(Float64, numSteps)

    for step in 1:numSteps
        U = disp[step]
        u_nodes = evaluate_at_grid_nodes(dh, U, :u)
        n_nodes = length(u_nodes)
        ux = [u_nodes[j][1] for j in 1:n_nodes]
        uy = [u_nodes[j][2] for j in 1:n_nodes]
        disp_points = [Point{2,Float64}(ux[j], uy[j]) for j in 1:n_nodes]
        UT[step] = disp_points
        UT_mag[step] = [sqrt(ux[j]^2 + uy[j]^2) for j in 1:n_nodes]
        ut_mag_max[step] = maximum(UT_mag[step])
    end
    return UT, UT_mag, ut_mag_max
end

UT, UT_mag, ut_mag_max = extract_solution(problem, grid)
numSteps = length(problem.post.times)


# Create displaced mesh per step
scale = 1.0

VV = [Point{2,Float64}(e[1], e[2]) for e in V]
VT = [VV .+ scale .* UT[i] for i in 1:numSteps]

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

# === Visualization setup ===
fig_disp = Figure(size=(1000,600))
stepStart = 1 # Start at undeformed
ax3 = Axis(fig_disp[1, 1], title = "Step: $stepStart")

xlims!(ax3, min_p[1], max_p[1])
ylims!(ax3, min_p[2], max_p[2])
hp = poly!(ax3, GeometryBasics.Mesh(VT[stepStart], F), 
               strokewidth = 0.5,
               color = UT_mag[stepStart], 
               transparency = false, 
               colormap = Reverse(:Spectral), 
               colorrange = (0, maximum(ut_mag_max)))

Colorbar(fig_disp[1, 2], hp.plots[1], label = "Displacement magnitude [mm]") 

incRange = 1:numSteps
hSlider = Slider(fig_disp[2, 1], range = incRange, startvalue = stepStart - 1, linewidth = 30)

on(hSlider.value) do stepIndex     
    hp[1] = GeometryBasics.Mesh(VT[stepIndex], F)
    hp.color = UT_mag[stepIndex]
    ax3.title = "Step: $stepIndex"
end

slidercontrol(hSlider, ax3)
display(GLMakie.Screen(), fig_disp)