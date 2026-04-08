using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics
using ComodoFerrite
using ComodoFerrite.Ferrite


## GLMakie setting 
GLMakie.closeall()

## Mesh 
boxDim = [10.0,40.0,10.0]
boxEl = [5, 10, 5]
E, V, F, Fb, Cb = hexbox(boxDim, boxEl)
grid = ComodoToFerrite(E, V)
Fb_back = Fb[Cb.==4]   
addface!(grid , "back", Fb_back) 

Fb_front = Fb[Cb.==3]  
addface!(grid , "front", Fb_front) 



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


function create_bc(dh, grid)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(grid, "back"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3]) # bcSupportList_Z
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



function assemble_element!(ke, ge, cell, cv, fv, mp, ue, ΓN, tn)
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
    for facet in 1:nfacets(cell)
        if (cellid(cell), facet) in ΓN
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
                T0 = Φ * tn
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
                        δT0 = δΦ * tn
                        ke[i, j] -= (δui ⋅ δT0) * dΓ0
                    end
                end
            end
        end
    end
end;

function assemble_global!(K, g, dh, cv, fv, mp, u, ΓN, tn)
    n = ndofs_per_cell(dh)
    ke = zeros(n, n)
    ge = zeros(n)

    # start_assemble resets K and g
    assembler = start_assemble(K, g)

    # Loop over all cells in the grid
    for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] # element dofs
        assemble_element!(ke, ge, cell, cv, fv, mp, ue, ΓN, tn)
        assemble!(assembler, global_dofs, ke, ge)
    end
    return
end;

function solve(E, ν, grid, traction_prescribed, numSteps)

    # --- Material ---
    μ = E / (2 * (1 + ν))
    λ = (E * ν) / ((1 + ν) * (1 - 2ν))
    mp = NeoHooke(μ, λ)

    # --- FEM setup ---
    dh = create_dofhandler(grid)
    dbcs = create_bc(dh,grid)
    cv, fv = create_values()

    ΓN = getfacetset(grid, "front")

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
    g = zeros(nd)

    # --- Parameters ---
    NEWTON_TOL = 1e-8
    NEWTON_MAXITER = 100

    Tf = 1.0                    # load factor
    Δt = Tf / numSteps

    # --- Initial condition ---
    apply!(un, dbcs)
    u .= un

    # --- Time stepping (UNCHANGED) ---
    for (step, t) in enumerate(0.0:Δt:Tf)
        println("\n=== Time step $step, t = $t ===")

        fill!(Δu, 0.0)
        newton_itr = 0

        # scaled traction vector
        tn = t .*traction_prescribed
        

        while true
            u .= un .+ Δu
            apply!(u, dbcs)

            assemble_global!(K, g, dh, cv, fv, mp, u, ΓN, tn)

            normg = norm(g[Ferrite.free_dofs(dbcs)])
            if normg < NEWTON_TOL
                println("  Converged in $newton_itr iterations")
                break
            elseif newton_itr ≥ NEWTON_MAXITER
                error("Newton failed to converge at time t = $t")
            end

            apply_zero!(K, g, dbcs)

            fill!(ΔΔu, 0.0)
            ΔΔu = K \ g
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


E = 10.0
ν = 0.3

traction_prescribed =  [0.0 , 0.0, -0.1]
numSteps = 10
UT, UT_mag, ut_mag_max = solve(E, ν, grid, traction_prescribed, numSteps)


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