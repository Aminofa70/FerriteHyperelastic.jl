using Revise
using FerriteHyperelastic
using Ferrite

input = InputStruct()
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Triangle, (nx, ny), corners)
    addfacetset!(grid, "support_1", x -> x[1] ≈ 0.0)
    addfacetset!(grid, "pressure", x -> x[1] ≈ Lx)
    return grid
end
# create fem values: interpolation and quadrature point
function create_values()
    ipu, ipp = Lagrange{RefTriangle,2}()^2, Lagrange{RefTriangle,1}() # P2/P1
    qr = QuadratureRule{RefTriangle}(3) # Higher quadrature order for stability
    qr_face = FacetQuadratureRule{RefTriangle}(1)
    return CellValues(qr, ipu), CellValues(qr, ipp), FacetValues(qr_face, ipu)
end
##################################################################
# create dof (degree of freedom)
function create_dofhandler(grid)
    ipu, ipp = Lagrange{RefTriangle,2}()^2, Lagrange{RefTriangle,1}()
    dh = DofHandler(grid)
    add!(dh, :u, ipu) # displacement dim = 3
    add!(dh, :p, ipp) # pressure dim = 1
    close!(dh)
    return dh
end;

# creat boundary condition Dirichlet (displacement)
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getfacetset(dh.grid, "support_1"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.close!(ch)
    return ch
end

####################################

function create_bc_force(dh)
    dbc = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(dbc, Ferrite.Dirichlet(:u, getfacetset(grid, "pressure"), (x, t) -> 0*x))
    Ferrite.close!(dbc)
    return dbc
end

function Ψ(F, p, μ, λ)
    Ic = tr(tdot(F))
    J = det(F)
    Js = (λ + p + sqrt((λ + p)^2.0 + 4.0 * λ * μ)) / (2.0 * λ)
    return p * (Js - J) + μ / 2 * (Ic - 3) - μ * log(Js) + λ / 2 * (Js - 1)^2
end;

# # Define the function to find Js numerically by solving f(Js) = 0
function constitutive_driver(F, p, μ, λ)
    # Compute all derivatives in one function call
    ∂²Ψ∂F², ∂Ψ∂F = Tensors.hessian(y -> Ψ(y, p, μ, λ), F, :all)
    ∂²Ψ∂p², ∂Ψ∂p = Tensors.hessian(y -> Ψ(F, y, μ, λ), p, :all)
    ∂²Ψ∂F∂p = Tensors.gradient(q -> Tensors.gradient(y -> Ψ(y, q, μ, λ), F), p)
    return ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂p, ∂²Ψ∂p², ∂²Ψ∂F∂p
end;

function make_constitutive_driver(μ, κ)
    return (F, p) -> constitutive_driver(F, p, μ, λ)
end

input.model_type = :plane_strain   # or :plane_strain; :plane_stress; :threeD
input.load_type = :traction

ν = 0.4999

μ = 80.1938
E = (2 * μ) * (1 + ν)
λ = (2 * μ * ν) / (1 - 2ν)
κ = E / (3 * (1 - 2 * ν))    # Bulk modulus

C10 = E / (4 * (1 + ν))
D1 = 6.0 * (1.0 - 2.0 * ν) / E

input.material = make_constitutive_driver(μ, λ)

# Define parameters for the plate and mesh
Lx, Ly = 3.17, 1.73  # Plate dimensions
nx, ny = 10, 10   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid

input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh)
# Create CellValues and FacetValues
input.cell_values_u, input.cell_values_p, input.facet_values = create_values()

input.ΓN = getfacetset(grid, "pressure")
input.facetsets = [input.ΓN]
input.traction = [20.2, 0.0]
input.tractions = Dict(1 => input.traction)
input.tol = 1e-6
## default
maxIterPerInc, totalTime, initInc, minInc, maxInc, totalInc = initialize_solver()
# change like the following if you need
# maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver(1000, 1.0, 1e-6,1e-15, 0.2,1000)

input.maxIterPerInc = maxIterPerInc
input.totalTime = totalTime
input.initInc = initInc
input.minInc = minInc
input.maxInc = maxInc
input.totalInc = totalInc

input.dof_F = []
input.dof_U = []

input.filename = "2D_Hyper"
input.output_dir= "/Users/aminalibakhshi/Desktop/vtu_geo/"
##################################################################
################  solution 

sol  = run_fem_hybrid(input)

dh = input.dh
U = sol.U_steps[end]
u_dof = dof_range(input.dh, :u)
p_dof = dof_range(dh, :p)


# Assume `dh` (DofHandler) and `U` (solution vector) are already defined
n_elements = getncells(input.dh.grid)  # Total number of elements

# Initialize global vectors for u and p
u_global = Float64[]
p_global = Float64[]

# Loop over all elements and collect u and p DOFs
for cell in CellIterator(dh)
    dofs = celldofs(cell)  # Get global DOF indices for this element
    u_local = U[dofs[u_dof]]
    p_local = U[dofs[p_dof]]
    # Append to global vectors
    append!(u_global, u_local)
    append!(p_global, p_local)
end

# Now `u_global` contains all displacement DOFs, and `p_global` all pressure DOFs
println("Total displacement DOFs: ", length(u_global))
println("Total pressure DOFs: ", length(p_global))

# Extract final displacement and evaluate at grid nodes
U_end = u_global
u_nodes = vec(evaluate_at_grid_nodes(input.dh, U_end, :u))
ux, uy = getindex.(u_nodes, 1), getindex.(u_nodes, 2)
@info "Max |ux| = $(maximum(abs, ux)), Max |uy| = $(maximum(abs, uy))"