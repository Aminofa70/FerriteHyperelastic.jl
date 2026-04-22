using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics
using ComodoFerrite
using ComodoFerrite.Ferrite
using FESolvers

const Vec = Ferrite.Vec
## GLMakie setting
GLMakie.closeall()

## Mesh
boxDim = [5, 5, 5]
boxEl = [10, 10, 10]
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

## Finite Element Values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefHexahedron, order}()^dim
    qr = QuadratureRule{RefHexahedron}(3)
    qr_face = FacetQuadratureRule{RefHexahedron}(2)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

## Create Degrees of freedom
function create_dofhandler(grid)
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefHexahedron, 1}()^3)
    close!(dh)
    return dh
end

function create_bc(dh)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "bottom"), (x, t) -> [0.0], [3]))
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "front"), (x, t) -> [0.0], [2]))
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "left"), (x, t) -> [0.0], [1]))
    close!(ch)
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
end

function assemble_element!(ke, ge, cell, cv, fv, mp, ue, ΓN, tn, time)
    reinit!(cv, cell)
    fill!(ke, 0.0)
    fill!(ge, 0.0)
    ndofs = getnbasefunctions(cv)

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        ∇u = function_gradient(cv, qp, ue)
        F = one(∇u) + ∇u
        C = tdot(F)
        S, ∂S∂C = constitutive_driver(C, mp)
        P = F ⋅ S
        I = one(S)
        ∂P∂F = otimesu(I, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(F', I)

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

    # Follower traction on Neumann boundary
    tn_current = time * tn
    for facet in 1:nfacets(cell)
        if FacetIndex(cellid(cell), facet) in ΓN
            reinit!(fv, cell, facet)
            for q_point in 1:getnquadpoints(fv)
                ∇u = function_gradient(fv, q_point, ue)
                F = one(∇u) + ∇u
                J = det(F)
                FinvT = inv(F)'
                C = F' ⋅ F
                Cinv = inv(C)
                N₀ = getnormal(fv, q_point)
                # α = N₀ ⋅ C⁻¹ ⋅ N₀
                m = Cinv ⋅ N₀
                α = N₀ ⋅ m
                sqrtα = sqrt(α)
                Φ = J * sqrtα
                T0 = Φ * tn_current
                dΓ0 = getdetJdV(fv, q_point)
                for i in 1:ndofs
                    δui = shape_value(fv, q_point, i)
                    ge[i] -= (δui ⋅ T0) * dΓ0
                    for j in 1:ndofs
                        ∇δuj = shape_gradient(fv, q_point, j)
                        δF = ∇δuj
                        δJ = J * (FinvT ⊡ δF)
                        δC = δF' ⋅ F + F' ⋅ δF
                        δα = -(m ⋅ (δC ⋅ m))
                        δΦ = sqrtα * δJ + (J / (2 * sqrtα)) * δα
                        δT0 = δΦ * tn_current
                        ke[i, j] -= (δui ⋅ δT0) * dΓ0
                    end
                end
            end
        end
    end
    return
end

function assemble_global!(K, g, dh, cv, fv, mp, u, ΓN, tn, time)
    n = ndofs_per_cell(dh)
    ke = zeros(n, n)
    ge = zeros(n)
    assembler = start_assemble(K, g)

    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs]
        assemble_element!(ke, ge, cell, cv, fv, mp, ue, ΓN, tn, time)
        assemble!(assembler, global_dofs, ke, ge)
    end
    return
end

## Post-processing storage
mutable struct NeoHookePost
    solutions::Vector{Vector{Float64}}
    times::Vector{Float64}
end

NeoHookePost() = NeoHookePost(Vector{Float64}[], Float64[])

## Problem struct
struct NeoHookeProblem{PD, PB, PP}
    def::PD
    buf::PB
    post::PP
end

struct NeoHookeModel{FS}
    dh::DofHandler
    ch::ConstraintHandler
    material::NeoHooke
    ΓN::FS
    tn::Vec{3, Float64}
end

struct NeoHookeBuffer{CV, FV, KT, T}
    cv::CV
    fv::FV
    K::KT
    r::Vector{T}
    u::Vector{T}
    time::Vector{T}
end

function build_buffer(model::NeoHookeModel)
    cv, fv = create_values()
    K = allocate_matrix(model.dh)
    r = zeros(ndofs(model.dh))
    u = zeros(ndofs(model.dh))
    return NeoHookeBuffer(cv, fv, K, r, u, [0.0])
end

function build_problem(def::NeoHookeModel)
    return NeoHookeProblem(def, build_buffer(def), NeoHookePost())
end

## FESolvers interface
FESolvers.getunknowns(p::NeoHookeProblem) = p.buf.u
FESolvers.getresidual(p::NeoHookeProblem) = p.buf.r
FESolvers.getjacobian(p::NeoHookeProblem) = p.buf.K

function FESolvers.update_to_next_step!(p::NeoHookeProblem, time)
    p.buf.time .= time
    return Ferrite.update!(p.def.ch, time)
end

function FESolvers.update_problem!(p::NeoHookeProblem, Δu, _)
    if !isnothing(Δu)
        apply_zero!(Δu, p.def.ch)
        p.buf.u .+= Δu
    end
    assemble_global!(
        p.buf.K, p.buf.r, p.def.dh, p.buf.cv, p.buf.fv,
        p.def.material, p.buf.u, p.def.ΓN, p.def.tn, p.buf.time[1]
    )
    return apply_zero!(p.buf.K, p.buf.r, p.def.ch)
end

function FESolvers.calculate_convergence_measure(p::NeoHookeProblem, args...)
    return norm(FESolvers.getresidual(p)[free_dofs(p.def.ch)])
end

function FESolvers.postprocess!(p::NeoHookeProblem, solver)
    push!(p.post.solutions, copy(p.buf.u))
    return push!(p.post.times, p.buf.time[1])
end

FESolvers.handle_converged!(::NeoHookeProblem) = nothing

## Material parameters
Emod = 1.0
ν = 0.3
μ = Emod / (2(1 + ν))
λ = Emod * ν / ((1 + ν) * (1 - 2ν))
material = NeoHooke(μ, λ)

## Setup
dh = create_dofhandler(grid)
ch = create_bc(dh)

## Neumann boundary (follower traction on top face)
ΓN = getfacetset(grid, "top")
traction = 0.6
tn = Vec{3}((0.0, 0.0, traction))

## Build and solve
def = NeoHookeModel(dh, ch, material, ΓN, tn)
problem = build_problem(def)

solver = QuasiStaticSolver(
    NewtonSolver(;
        linsolver = BackslashSolver(), linesearch = NoLineSearch(),
        maxiter = 30, tolerance = 1.0e-9,
        update_jac_first = true, update_jac_each = true
    )
    ,
    FixedTimeStepper(; num_steps = 10, Δt = 0.1, t_start = 0.0)
)

solve_problem!(problem, solver)

function solution(disp, numSteps)
    UT = Vector{Vector{Point{3, Float64}}}(undef, numSteps)
    UT_mag = Vector{Vector{Float64}}(undef, numSteps)
    ut_mag_max = zeros(Float64, numSteps)

    dh = create_dofhandler(grid)

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

disp = problem.post.solutions
numSteps = length(problem.post.times)

UT, UT_mag, ut_mag_max = solution(disp, numSteps)

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
