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

Vec = Ferrite.Vec
## GLMakie setting
GLMakie.closeall()

## Mesh
sampleSize = 10.0
pointSpacing = 2.0
strainApplied = 0.5
loadingOption = "tension"

boxDim = sampleSize .* [1, 1, 1] #Dimensionsions for the box in each direction
boxEl = ceil.(Int64,boxDim./pointSpacing) # Number of elements to use in each direction
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

if loadingOption == "tension"
    displacement_prescribed = strainApplied * sampleSize
elseif loadingOption == "compression"
    displacement_prescribed = -strainApplied * sampleSize
end

## Finite Element Values
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
    dbc = Dirichlet(:u, getfacetset(dh.grid, "bottom"), (x, t) -> [0.0], [3])
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(dh.grid, "front"), (x, t) -> [0.0], [2])
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(dh.grid, "left"), (x, t) -> [0.0], [1])
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(dh.grid, "top"), (x, t) -> [displacement_prescribed * t], [3])
    add!(ch, dbc)
    Ferrite.close!(ch)
    Ferrite.update!(ch, 0.0)
    return ch
end

## Material definition: Neo-Hookean + Exponential Fiber
struct NeoHookeFiber
    μ::Float64
    λ::Float64
    ξ::Float64
    α::Float64
    β::Float64
    θ::Float64  # in RADIANS
    ϕ::Float64  # in RADIANS
end

# Matrix-only strain energy (Neo-Hookean)
function Ψ_matrix(C, mp::NeoHookeFiber)
    μ = mp.μ
    λ = mp.λ
    Ic = tr(C)
    J = sqrt(det(C))
    return μ / 2 * (Ic - 3) - μ * log(J) + λ / 2 * (log(J))^2
end

# Fiber strain energy only
# Handles α = 0 limit: Ψ_fib = (ξ/β)(Iₙ - 1)^β
# Handles α > 0:        Ψ_fib = (ξ/(αβ))(exp[α(Iₙ - 1)^β] - 1)
function Ψ_fiber(C, mp::NeoHookeFiber)
    n_r = Vec{3}((cos(mp.θ) * sin(mp.ϕ), sin(mp.θ) * sin(mp.ϕ), cos(mp.ϕ)))
    Iₙ = n_r ⋅ (C ⋅ n_r)
    λ₀ = 1.0
    I₀= sqrt(λ₀)
    x = max(Iₙ - I₀, 0.0)  # tension only (Heaviside built in)
    if mp.α ≈ 0.0
        # Power law limit: lim α→0 of (ξ/(αβ))(exp(α x^β) - 1) = (ξ/β) x^β
        return (mp.ξ / mp.β) * x^mp.β
    else
        return (mp.ξ / (mp.α * mp.β)) * (exp(mp.α * x^mp.β) - 1.0)
    end
end

# Full strain energy
function Ψ(C, mp::NeoHookeFiber)
    return Ψ_matrix(C, mp) + Ψ_fiber(C, mp)
end

# Constitutive driver
function constitutive_driver(C, mp::NeoHookeFiber)
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, mp), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end

function assemble_element!(ke, ge, cell, cv, mp, ue)
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
end

function assemble_global!(K, g, dh, cv, mp, u)
    n = ndofs_per_cell(dh)
    ke = zeros(n, n)
    ge = zeros(n)
    assembler = start_assemble(K, g)
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs]
        assemble_element!(ke, ge, cell, cv, mp, ue)
        assemble!(assembler, global_dofs, ke, ge)
    end
    return
end

## Material parameters
E_mod = 1.0
ν = 0.4
μ = E_mod / (2 * (1 + ν))
λ = (E_mod * ν) / ((1 + ν) * (1 - 2ν))

# Fiber parameters
ξ = 1.0     # fiber stiffness scaling (ξ > 0)
α = 2.0     # exponential coefficient (α ≥ 0; α=0 gives power law)
β = 2.0     # power exponent (β ≥ 2; β > 2 for smooth transition)
θ = deg2rad(0.0)    # FEBio: <theta>0</theta>
ϕ = deg2rad(45.0)   # FEBio: <phi>90</phi>


mp = NeoHookeFiber(μ, λ, ξ, α, β, θ, ϕ)

## Post-processing storage
mutable struct NeoHookePost
    solutions::Vector{Vector{Float64}}
    times::Vector{Float64}
end
NeoHookePost() = NeoHookePost(Vector{Float64}[], Float64[])

## Problem struct
struct NeoHookeProblem{PD,PB,PP}
    def::PD
    buf::PB
    post::PP
end

struct NeoHookeModel{DH,CH,M}
    dh::DH
    ch::CH
    material::M
end

function NeoHookeModel()
    dh = create_dofhandler(grid)
    ch = create_bc(dh)
    return NeoHookeModel(dh, ch, mp)
end

struct NeoHookeBuffer{CV,KT,T}
    cv::CV
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
    return NeoHookeBuffer(cv, K, r, u, [0.0])
end

build_problem(def::NeoHookeModel) = NeoHookeProblem(def, build_buffer(def), NeoHookePost())

## FESolvers interface
FESolvers.getunknowns(p::NeoHookeProblem) = p.buf.u
FESolvers.getresidual(p::NeoHookeProblem) = p.buf.r
FESolvers.getjacobian(p::NeoHookeProblem) = p.buf.K

function FESolvers.update_to_next_step!(p::NeoHookeProblem, time)
    p.buf.time .= time
    Ferrite.update!(p.def.ch, time)
    apply!(p.buf.u, p.def.ch)
end

function FESolvers.update_problem!(p::NeoHookeProblem, Δu, _)
    if !isnothing(Δu)
        apply_zero!(Δu, p.def.ch)
        p.buf.u .+= Δu
    end
    assemble_global!(p.buf.K, p.buf.r, p.def.dh, p.buf.cv, p.def.material, p.buf.u)
    apply_zero!(p.buf.K, p.buf.r, p.def.ch)
end

FESolvers.calculate_convergence_measure(p::NeoHookeProblem, args...) = norm(FESolvers.getresidual(p)[free_dofs(p.def.ch)])

function FESolvers.postprocess!(p::NeoHookeProblem, solver)
    push!(p.post.solutions, copy(p.buf.u))
    push!(p.post.times, p.buf.time[1])
    println("Time step $(length(p.post.times)) completed, t = $(p.buf.time[1])")
end

FESolvers.handle_converged!(::NeoHookeProblem) = nothing

## Build and solve
def = NeoHookeModel()
problem = build_problem(def)

solver = QuasiStaticSolver(
    NewtonSolver(;
        linsolver=BackslashSolver(),
        linesearch=NoLineSearch(),
        maxiter=30,
        tolerance=1.0e-9,
        update_jac_first=true,
        update_jac_each=true),
    FixedTimeStepper(; num_steps=10, Δt=0.1, t_start=0.0)
)
solve_problem!(problem, solver)

## Post-processing
function solution(disp, numSteps)
    UT = Vector{Vector{Point{3,Float64}}}(undef, numSteps)
    UT_mag = Vector{Vector{Float64}}(undef, numSteps)
    ut_mag_max = zeros(Float64, numSteps)

    dh = create_dofhandler(grid)

    for step in 1:numSteps
        U = disp[step]
        u_nodes = vec(evaluate_at_grid_nodes(dh, U, :u))
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

disp = problem.post.solutions
numSteps = length(problem.post.times)

UT, UT_mag, ut_mag_max = solution(disp, numSteps)

numInc = length(UT)

# Create displaced mesh per step
scale = 1.0
VT = [V .+ scale .* UT[i] for i in 1:(numSteps)]
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