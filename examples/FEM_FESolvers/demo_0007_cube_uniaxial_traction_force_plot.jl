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

GLMakie.closeall()

## Mesh
sampleSize = 10.0
pointSpacing = 2.0
strainApplied = 0.5 # Equivalent linear strain
loadingOption = "compression" # "tension" or "compression"
# Creating a hexahedral mesh for a cube
boxDim = sampleSize .* [1, 1, 1] # Dimensionsions for the box in each direction
boxEl = ceil.(Int64, boxDim ./ pointSpacing) # Number of elements to use in each direction
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
    return
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

#########################################################
## Internal force via Bᵀσ (general formulation)
#########################################################
function assemble_internal_forces(dh, cv, mp, u)
    f_int = zeros(ndofs(dh))
    n = getnbasefunctions(cv)
    fe = zeros(n)

    for cell in CellIterator(dh)
        reinit!(cv, cell)
        fill!(fe, 0.0)
        ue = u[celldofs(cell)]

        for qp in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, qp)
            ∇u = function_gradient(cv, qp, ue)
            F = one(∇u) + ∇u
            C = tdot(F)
            S, _ = constitutive_driver(C, mp)
            P = F ⋅ S

            for i in 1:n
                ∇Ni = shape_gradient(cv, qp, i)
                fe[i] += (∇Ni ⊡ P) * dΩ
            end
        end

        assemble!(f_int, celldofs(cell), fe)
    end
    return f_int
end

#########################################################
## Boundary nodes and reaction extraction
#########################################################
function get_boundary_nodes(grid, facetset)
    boundary_nodes = Set{Int}()
    for (cell_id, facet_id) in facetset
        cell = getcells(grid)[cell_id]
        facet_node_ids = Ferrite.facets(cell)[facet_id]
        for n in facet_node_ids
            push!(boundary_nodes, n)
        end
    end
    return boundary_nodes
end

function compute_reaction(f_int, boundary_nodes)
    Rx = Float64[]
    Ry = Float64[]
    Rz = Float64[]
    for node in boundary_nodes
        push!(Rx, f_int[3 * (node - 1) + 1])
        push!(Ry, f_int[3 * (node - 1) + 2])
        push!(Rz, f_int[3 * (node - 1) + 3])
    end
    return Rx, Ry, Rz
end

## Material
E_mod = 1.0
ν = 0.4
μ = E_mod / (2 * (1 + ν))
λ = (E_mod * ν) / ((1 + ν) * (1 - 2ν))
mp = NeoHooke(μ, λ)

## Post-processing storage
mutable struct NeoHookePost
    solutions::Vector{Vector{Float64}}
    times::Vector{Float64}
    reactions::Vector{Vector{Float64}}
end
NeoHookePost() = NeoHookePost(Vector{Float64}[], Float64[], Vector{Float64}[])

## Problem struct
struct NeoHookeProblem{PD, PB, PP}
    def::PD
    buf::PB
    post::PP
end

struct NeoHookeModel{DH, CH, M}
    dh::DH
    ch::CH
    material::M
end

function NeoHookeModel()
    dh = create_dofhandler(grid)
    ch = create_bc(dh)
    return NeoHookeModel(dh, ch, mp)
end

struct NeoHookeBuffer{CV, KT, T}
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
    return apply!(p.buf.u, p.def.ch)
end

function FESolvers.update_problem!(p::NeoHookeProblem, Δu, _)
    if !isnothing(Δu)
        apply_zero!(Δu, p.def.ch)
        p.buf.u .+= Δu
    end
    assemble_global!(p.buf.K, p.buf.r, p.def.dh, p.buf.cv, p.def.material, p.buf.u)
    return apply_zero!(p.buf.K, p.buf.r, p.def.ch)
end

FESolvers.calculate_convergence_measure(p::NeoHookeProblem, args...) = norm(FESolvers.getresidual(p)[free_dofs(p.def.ch)])

function FESolvers.postprocess!(p::NeoHookeProblem, solver)
    push!(p.post.solutions, copy(p.buf.u))
    push!(p.post.times, p.buf.time[1])

    # Compute f_int (general formulation)
    f_int = assemble_internal_forces(p.def.dh, p.buf.cv, p.def.material, p.buf.u)
    return push!(p.post.reactions, copy(f_int))
end

FESolvers.handle_converged!(::NeoHookeProblem) = nothing

## Build and solve
def = NeoHookeModel()
problem = build_problem(def)

solver = QuasiStaticSolver(
    NewtonSolver(; tolerance = 1.0e-6, maxiter = 30),
    FixedTimeStepper(; num_steps = 10, Δt = 0.1, t_start = 0.0)
)

solve_problem!(problem, solver)

#########################################################
## Post-processing: displacements
#########################################################
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

#########################################################
## Reaction force curves
#########################################################
top_nodes = get_boundary_nodes(grid, getfacetset(grid, "top"))

Fz_curve = zeros(numSteps)
time_curve = zeros(numSteps)

for i in 1:numSteps
    Rx, Ry, Rz = compute_reaction(problem.post.reactions[i], top_nodes)
    Fz_curve[i] = sum(Rz)
    time_curve[i] = problem.post.times[i]
end

#########################################################
## Visualization (FEBio.jl style)
#########################################################
numInc = length(UT)   # = numSteps + 1

scale = 1.0
VT = [V .+ scale .* UT[i] for i in 1:numInc]   #

min_p = minp([minp(Vi) for Vi in VT])
max_p = maxp([maxp(Vi) for Vi in VT])

incRange = 0:(numInc - 1)   # 0 → 10

fig = Figure(size = (1600, 800))
stepStart = 1

ax1 = AxisGeom(
    fig[1, 1], title = "Step: $stepStart",
    limits = (min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3])
)

hp1 = meshplot!(
    ax1, Fb, VT[stepStart + 1];   # shift index
    strokewidth = 2,
    color = UT_mag[stepStart + 1],
    transparency = false,
    colormap = Reverse(:Spectral),
    colorrange = (0.0, maximum(ut_mag_max))
)

Colorbar(fig[1, 2], hp1.plots[1], label = "Displacement magnitude [mm]")

ax3 = Axis(
    fig[1, 3], title = "Step: $stepStart", aspect = AxisAspect(1),
    xlabel = "Time [s]", ylabel = "Force [N]"
)

lines!(ax3, time_curve, Fz_curve, color = :red, linewidth = 3)

hp3 = scatter!(
    ax3,
    Point{2, Float64}(time_curve[stepStart + 1], Fz_curve[stepStart + 1]),
    markersize = 15, color = :red
)

ax4 = Axis(
    fig[1, 4], title = "Step: $stepStart", aspect = AxisAspect(1),
    xlabel = "Max displacement [mm]", ylabel = "Force [N]"
)

lines!(ax4, ut_mag_max, Fz_curve, color = :blue, linewidth = 3)

hp4 = scatter!(
    ax4,
    Point{2, Float64}(ut_mag_max[stepStart + 1], Fz_curve[stepStart + 1]),
    markersize = 15, color = :blue
)

hSlider = Slider(fig[2, :], range = incRange, startvalue = stepStart, linewidth = 30)

on(hSlider.value) do stepIndex
    i = stepIndex + 1   # convert 0-based → 1-based

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
