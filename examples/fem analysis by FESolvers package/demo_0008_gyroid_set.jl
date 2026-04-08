using FESolvers
using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics
using IterativeSolvers
using ComodoFerrite
using ComodoFerrite.Ferrite
using SparseArrays , LinearAlgebra
using Geogram

GLMakie.closeall()

###### 
# Control parameters 
E_youngs = 1.0
ν = 0.3

displacement_prescribed = 0.5

###### 
pointSpacing = 0.5
level = 0.75
s1 = -0.1*π
s2 =  4.1*π
ns = spacing2numsteps(s2-s1, pointSpacing/3.0; close_loop=false)
x = range(s1, s2, ns)

type=:G
F1, V1 = tpms(type; x=x, level= -level, cap = true, padValue=1e8, side=:positive)
F2, V2 = tpms(type; x=x, level= level, cap = true, padValue=1e8, side=:negative)

np1 = spacing2numvertices(F1, V1, pointSpacing)
F1,V1 = ggremesh(F1,V1; nb_pts=np1)

np2 = spacing2numvertices(F2, V2, pointSpacing)
F2,V2 = ggremesh(F2,V2; nb_pts=np2)

V1_in = faceinteriorpoint(F1,V1,1)
V2_in = faceinteriorpoint(F2,V2,1)

Fb, Vb, Cb = joingeom(F1, V1, F2, V2)

V_regions = [V1_in, V2_in]

vol1 = pointSpacing^3 / (6.0*sqrt(2.0))
region_vol = [vol1, vol1]
stringOpt = "paAqYQ"
E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb; facetmarkerlist=Cb, V_regions=V_regions, region_vol=region_vol, stringOpt)

grid = ComodoToFerrite(E, V)
# Find boundary condition faces 
indicesTopNodes = Vector{Int}()
indicesBottomNodes = Vector{Int}()
Z = [v[3] for v in V]
zMax = maximum(Z)
zMin = minimum(Z)
search_tol = pointSpacing/2.0
for f in Fb    
    bTop = Z[f].>(zMax-search_tol)
    bBottom = Z[f].<(zMin+search_tol)
    if all(bTop)
        append!(indicesTopNodes, f)        
    end
    if all(bBottom)
        append!(indicesBottomNodes, f)        
    end
end
indicesTopNodes = unique(indicesTopNodes)
indicesBottomNodes = unique(indicesBottomNodes)

addfacetset!(grid, "top", boundary_facets(grid, indicesTopNodes))
addfacetset!(grid, "bottom", boundary_facets(grid, indicesBottomNodes))


## Finite Elemenet Values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefTetrahedron,order}()^dim
    qr = QuadratureRule{RefTetrahedron}(2)
    qr_face = FacetQuadratureRule{RefTetrahedron}(1)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

## Create Degrees of freedom
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefTetrahedron,1}()^3)
    Ferrite.close!(dh)
    return dh
end

function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(dh.grid, "bottom"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3])
    add!(ch, dbc)
    # Scale: t goes 0→1, displacement goes 0→max_displacement
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


μ = E_youngs / (2 * (1 + ν))
λ = (E_youngs * ν) / ((1 + ν) * (1 - 2ν))
mp = NeoHooke(μ, λ)

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
    println("Step $(length(p.post.times)): t = $(p.buf.time[1]), max|u| = $(maximum(abs, p.buf.u))")
end

FESolvers.handle_converged!(::NeoHookeProblem) = nothing

## Build and solve
def = NeoHookeModel()
problem = build_problem(def)


num_steps = 10
Δt = 1/num_steps

solver = QuasiStaticSolver(
    NewtonSolver(;
        linsolver=BackslashSolver(),
        linesearch=NoLineSearch(),
        maxiter= 100,
        tolerance= 1.0e-6,
        update_jac_first=true,
        update_jac_each=true),
    FixedTimeStepper(; num_steps=num_steps, Δt= Δt, t_start=0.0)
)

solve_problem!(problem, solver)



function solution(disp , numSteps)
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

incRange = 1:numSteps
# Create displaced mesh per step
scale = 1.0
VT = [V .+ scale .* UT[i] for i in 1:numSteps]

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

fig2 = Figure(size=(1200, 800))
stepStart = incRange[end]
ax1 = AxisGeom(fig2[1, 1], title = "Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp1 = meshplot!(ax1, Fb, VT[end]; strokewidth=1, color=UT_mag[end], transparency=false, colormap = Reverse(:Spectral), colorrange=(0.0, maximum(ut_mag_max)))
Colorbar(fig2[1, 2], hp1.plots[1], label = "Displacement magnitude [mm]") 


hSlider2 = Slider(fig2[2, :], range = incRange, startvalue = stepStart,linewidth=30)

on(hSlider2.value) do stepIndex 
    hp1[1] = GeometryBasics.Mesh(VT[stepIndex], Fb)
    hp1.color = UT_mag[stepIndex]
    ax1.title = "Step: $stepIndex"
end

screen = display(GLMakie.Screen(), fig2)
GLMakie.set_title!(screen, "FEBio example")