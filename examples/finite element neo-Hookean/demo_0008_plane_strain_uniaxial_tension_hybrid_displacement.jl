using Revise
using FerriteHyperelastic
using Ferrite

input = InputStruct()
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "support_1", x -> x[1] ≈ 0.0)
    addfacetset!(grid, "disp", x -> x[1] ≈ Lx)
    return grid
end
function create_values()
    ipu, ipp = Lagrange{RefQuadrilateral,2}()^2, Lagrange{RefQuadrilateral,1}() # P2/P1
    qr = QuadratureRule{RefQuadrilateral}(3) # Higher quadrature order for stability
    qr_face = FacetQuadratureRule{RefQuadrilateral}(1)
    return CellValues(qr, ipu), CellValues(qr, ipp), FacetValues(qr_face, ipu)
end

function create_dofhandler(grid)
    ipu, ipp = Lagrange{RefQuadrilateral,2}()^2, Lagrange{RefQuadrilateral,1}()
    dh = DofHandler(grid)
    add!(dh, :u, ipu) # displacement dim = 3
    add!(dh, :p, ipp) # pressure dim = 1
    close!(dh)
    return dh
end;

function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "support_1"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getfacetset(dh.grid, "disp"), (x, t) -> [t], [1]))
    Ferrite.close!(ch)
    Ferrite.update!(ch, 0.0)
    return ch
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
input.load_type = :displacement

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
nx, ny = 4, 4 # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid

input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh)
# Create CellValues and FacetValues
input.cell_values_u, input.cell_values_p, input.facet_values = create_values()

##################################

input.dof_F = []
input.dof_U = []


input.displacement = 0.5
input.tol = 1e-6

input.filename = "2D_Hyper"
input.output_dir= "/Users/aminalibakhshi/Desktop/vtu_geo/"

#maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver()

maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver(500,1.0,1e-3,1e-15,0.1,1000)
input.maxIterPerInc = maxIterPerInc
input.totalTime = totalTime
input.initInc = initInc
input.minInc = minInc
input.maxInc = maxInc
input.totalInc = totalInc

sol = run_fem_hybrid(input)
# Final displacement vector
U = sol.U_steps[end]

# Displacement at all nodes as a Vector{Vec{2,Float64}}
u_nodes = vec(evaluate_at_grid_nodes(input.dh, U, :u))

# Split into scalar components
ux = getindex.(u_nodes, 1)
uy = getindex.(u_nodes, 2)
@info "Max |ux| = $(maximum(abs.(ux)))"
@info "Max |uy| = $(maximum(abs.(uy)))"