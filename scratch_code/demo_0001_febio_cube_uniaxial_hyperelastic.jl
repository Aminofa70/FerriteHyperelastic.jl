using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics
using ComodoFerrite
using ComodoFerrite.Ferrite
using TimerOutputs, ProgressMeter, IterativeSolvers
##########################################################################################
## Mesh 
GLMakie.closeall()

boxDim = [10, 10, 10]
boxEl = [5, 5, 5]
E, V, F, Fb, Cb = hexbox(boxDim, boxEl)
grid = ComodoToFerrite(E, V, Ferrite.Hexahedron; Fb, Cb)
##########################################################################################
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
##########################################################################################
## Dof handler 
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefHexahedron,1}()^3)
    Ferrite.close!(dh)
    return dh
end
##########################################################################################
function create_bc(dh, grid)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> [0.0], [3]) # bcSupportList_Z
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(grid, "front"), (x, t) -> [0.0], [2]) # bcSupportList_Y
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0], [1]) # bcSupportList_X
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(grid, "top"), (x, t) -> [t], [3]) # bcPrescribeList_Z
    add!(ch, dbc)
    Ferrite.close!(ch)
    Ferrite.update!(ch, 0.0)
    return ch
end
##########################################################################################
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
    # Compute all derivatives in one function call
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, mp), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end;

function assemble_element!(ke, ge, cell, cv, mp, ue)
    # Reinitialize cell values, and reset output arrays
    reinit!(cv, cell)
    fill!(ke, 0.0)
    fill!(ge, 0.0)

    ndofs = getnbasefunctions(cv)

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        # Compute deformation gradient F and right Cauchy-Green tensor C
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
            # Test function and gradient
            ∇δui = shape_gradient(cv, qp, i)
            # Add contribution to the residual from this test function
            ge[i] += (∇δui ⊡ P) * dΩ

            ∇δui∂P∂F = ∇δui ⊡ ∂P∂F # Hoisted computation
            for j in 1:ndofs
                ∇δuj = shape_gradient(cv, qp, j)
                # Add contribution to the tangent
                ke[i, j] += (∇δui∂P∂F ⊡ ∇δuj) * dΩ
            end
        end
    end
    return
end;

function assemble_global!(K, g, dh, cv, mp, u)
    n = ndofs_per_cell(dh)
    ke = zeros(n, n)
    ge = zeros(n)

    # start_assemble resets K and g
    assembler = start_assemble(K, g)

    # Loop over all cells in the grid
    @timeit "assemble" for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] # element dofs
        @timeit "element assemble" assemble_element!(ke, ge, cell, cv, mp, ue)
        assemble!(assembler, global_dofs, ke, ge)
    end
    return
end;

dh = create_dofhandler(grid)
ch = create_bc(dh, grid)
# Create CellValues and FacetValues
cv, _ = create_values()

# Material parameters
E = 1.0
ν = 0.4
μ = E / (2(1 + ν))
λ = (E * ν) / ((1 + ν) * (1 - 2ν))
mp = NeoHooke(μ, λ)

dbcs = ch
# Pre-allocation of vectors for the solution and Newton increments
_ndofs = ndofs(dh)
un = zeros(_ndofs) # previous solution vector
u = zeros(_ndofs)
Δu = zeros(_ndofs)
ΔΔu = zeros(_ndofs)
apply!(un, dbcs)

# Create sparse matrix and residual vector
K = allocate_matrix(dh)
g = zeros(_ndofs)

# Perform Newton iterations

NEWTON_TOL = 1.0e-8
NEWTON_MAXITER = 30

Tf = 10.0
Δt = 0.1

# Initialize solution
un = zeros(ndofs(dh))  # Previous converged solution
u = zeros(ndofs(dh))   # Current solution
Δu = zeros(ndofs(dh))  # Current increment
ΔΔu = zeros(ndofs(dh)) # Newton correction

# Apply initial boundary conditions (at t=0)
apply!(un, dbcs)
u .= un  # Start from initial condition

for (step, t) in enumerate(0.0:Δt:Tf)
    println("\n=== Time step $step, t = $t ===")

    # Update boundary conditions for current time
    Ferrite.update!(dbcs, t)

    # Start Newton iteration for this time step
    newton_itr = 0
    fill!(Δu, 0.0)  # Reset increment for this time step
    converged = false

    # Newton-Raphson iteration
    while !converged && newton_itr < NEWTON_MAXITER
        # Current displacement guess = previous converged + current increment
        u .= un .+ Δu

        # Ensure Dirichlet BCs are satisfied exactly
        apply!(u, dbcs)
        apply!(Δu, dbcs)

        # # Compute residual and tangent
        # K = allocate_matrix(dh)
        # g = zeros(ndofs(dh))
        assemble_global!(K, g, dh, cv, mp, u)

        # Apply boundary conditions to system
        apply_zero!(K, g, dbcs)

        # Check convergence
        normg = norm(g)
        println("  Newton iter $newton_itr, |g| = $normg")

        if normg < NEWTON_TOL
            converged = true
            println("  Converged in $newton_itr iterations")
            break
        end

        # Solve for Newton correction
        fill!(ΔΔu, 0.0)
        IterativeSolvers.cg!(ΔΔu, K, g; maxiter=1000)
        apply_zero!(ΔΔu, dbcs)

        # Update displacement increment
        Δu .-= ΔΔu
        newton_itr += 1
    end

    if !converged
        error("Newton failed to converge at time t = $t")
    end

    # Update converged solution for next time step
    un .= u

    # Optional: Print displacement info
    println("  Max displacement: $(maximum(abs.(u)))")
end

# Save ONLY the last time step
println("\nSaving final solution at t = $Tf")
@timeit "export" begin
    VTKGridFile("hyperelasticity_final", dh) do vtk
        write_solution(vtk, dh, u)
    end
end

println("\n=== Analysis complete ===")
print_timer(title="Analysis with $(getncells(grid)) elements", linechars=:ascii)
return u

maximum(u)
minimum(u)

