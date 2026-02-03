using NonlinearSolve
using LinearAlgebra
using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics
using IterativeSolvers
using ComodoFerrite
using ComodoFerrite.Ferrite
using LinearSolve

GLMakie.closeall()

pointSpacing = 2.0
strainApplied = 0.5 # Equivalent linear strain
loadingOption ="tension" # "tension" or "compression"

E_youngs = 1
ν =0.4

###### 
# Creating a hexahedral mesh for a cube 
boxDim = [5.0, 20.0, 60.0] # Dimensionsions for the box in each direction
boxEl = ceil.(Int64,boxDim./pointSpacing) # Number of elements to use in each direction 
E,V,F,Fb,Cb = hexbox(boxDim,boxEl)

grid = ComodoToFerrite(E, V)
Fb_bottom = Fb[Cb.==1]
addface!(grid , "bottom", Fb_bottom) 

Fb_front = Fb[Cb.==3]  
addface!(grid , "front", Fb_front) 

Fb_top = Fb[Cb.==2] 
addface!(grid , "top", Fb_top)   

Fb_left = Fb[Cb.==6]
addface!(grid , "left", Fb_left)   

# Defining displacement of the top surface in terms of x, y, and z components
if loadingOption=="tension"
    displacement_prescribed = strainApplied*boxDim[3]
elseif loadingOption=="compression"
    displacement_prescribed = -strainApplied*boxDim[3]
end

## Finite Elemenet Values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefHexahedron,order}()^dim
    qr = QuadratureRule{RefHexahedron}(2)
    qr_face = FacetQuadratureRule{RefHexahedron}(1)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

## Create Degrees of freedom
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefHexahedron,1}()^3)
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


function assemble_global_Jacobian!(K, dh, cv, mp, u)
    n = ndofs_per_cell(dh)
    ke = zeros(n, n)

    assembler = start_assemble(K)

    # Loop over all cells in the grid
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] # element dofs
        Ferrite.reinit!(cv, cell)
        fill!(ke, 0.0)
        ndofs = getnbasefunctions(cv)
        for qp in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, qp)
            ∇u = function_gradient(cv, qp, ue)
            F = one(∇u) + ∇u
            C = tdot(F) # F' ⋅ F
            S, ∂S∂C = constitutive_driver(C, mp)
            I = one(S)
            ∂P∂F = otimesu(I, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(F', I)

            # Loop over test functions
            for i in 1:ndofs
                ∇δui = shape_gradient(cv, qp, i)
                ∇δui∂P∂F = ∇δui ⊡ ∂P∂F
                for j in 1:ndofs
                    ∇δuj = shape_gradient(cv, qp, j)
                    ke[i, j] += (∇δui∂P∂F ⊡ ∇δuj) * dΩ
                end
            end
        end
        assemble!(assembler, global_dofs, ke)
    end
    return K
end;

function assemble_global_residual!(g, dh, cv, mp, u)
    fill!(g, 0.0)
    n = ndofs_per_cell(dh)
    ge = zeros(n)
    # Loop over all cells in the grid
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] # element dofs
        Ferrite.reinit!(cv, cell)
        fill!(ge, 0.0)
        ndofs = getnbasefunctions(cv)
        for qp in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, qp)
            ∇u = function_gradient(cv, qp, ue)
            F = one(∇u) + ∇u
            C = tdot(F) # F' ⋅ F
            S, _ = constitutive_driver(C, mp)
            P = F ⋅ S
            # Loop over test functions
            for i in 1:ndofs
                ∇δui = shape_gradient(cv, qp, i)
                ge[i] += (∇δui ⊡ P) * dΩ
            end
        end
        assemble!(g, global_dofs, ge)
    end
    return g
end;

function f_nonlinear!(res, u, p)
    dh, mp, dbcs = p
    cv, _ = create_values()

    assemble_global_residual!(res, dh, cv, mp, u)
    apply_zero!(res, dbcs)
end


function j_nonlinear!(K, u, p)
    dh, mp, dbcs = p
    cv, _ = create_values()

    fill!(K.nzval, 0.0)
    assemble_global_Jacobian!(K, dh, cv, mp, u)
    apply!(K, dbcs)
end

function solve_nonlinear(E, ν, grid, displacement_prescribed, numSteps)

    # --- Material ---
    μ = E / (2 * (1 + ν))
    λ = (E * ν) / ((1 + ν) * (1 - 2 * ν))   
    mp = NeoHooke(μ, λ)

    # --- FEM setup ---
    dh = create_dofhandler(grid)
    dbcs = create_bc(dh)

    UT = Vector{Vector{Point{3,Float64}}}(undef, numSteps + 1)
    UT_mag = Vector{Vector{Float64}}(undef, numSteps + 1)
    ut_mag_max = zeros(Float64, numSteps + 1)

    # --- State ---
    u = zeros(ndofs(dh))

    params = (dh, mp, dbcs)

    # Jacobian structure
    K_proto = allocate_matrix(dh)

    Tf = displacement_prescribed
    Δt = Tf / numSteps

    for (step, t) in enumerate(0.0:Δt:Tf)
        println("\n=== Time step $step, t = $t ===")

        Ferrite.update!(dbcs, t)
        apply!(u, dbcs)

        nl_func = NonlinearFunction(
            f_nonlinear!;
            jac = j_nonlinear!,
            jac_prototype = K_proto
        )

        prob = NonlinearProblem(nl_func, u, params)

        sol = NonlinearSolve.solve(
            prob,
            NewtonRaphson(linesearch = BackTracking());
            abstol = 1e-8,
            reltol = 1e-8,
            maxiters = 30
        )

        @assert sol.retcode == SciMLBase.ReturnCode.Success
        u .= sol.u   

        # --- Postprocessing ---
        u_nodes = vec(evaluate_at_grid_nodes(dh, u, :u))
        disp_points = [Point{3,Float64}(u_nodes[j]) for j in eachindex(u_nodes)]

        UT[step] = disp_points
        UT_mag[step] = norm.(disp_points)
        ut_mag_max[step] = maximum(UT_mag[step])
    end

    return UT, UT_mag, ut_mag_max
end


numSteps = 10
UT, UT_mag, ut_mag_max = solve_nonlinear(E_youngs, ν, grid, displacement_prescribed, numSteps)


# Create displaced mesh per step
#scale = 1.0
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
               shading = true,
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