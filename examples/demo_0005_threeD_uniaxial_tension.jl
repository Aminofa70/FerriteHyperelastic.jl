using Revise
using FerriteHyperelastic
using Ferrite
using PUPM
input = InputStruct()
# Function to create a 3D grid
function create_grid(Lx, Ly, Lz, nx, ny, nz)
    left = Ferrite.Vec(0.0, 0.0, 0.0)
    right = Ferrite.Vec(Lx, Ly, Lz)
    grid = generate_grid(Hexahedron, (nx, ny, nz), left, right)
    return grid
end
# Function to create cell and facet values
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
# Function to create a DOF handler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefHexahedron, 1}()^3)
    Ferrite.close!(dh)
    return dh
end
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getfacetset(dh.grid, "left"), (x, t) -> [0.0, 0.0, 0.0], [1, 2,3]))
    Ferrite.close!(ch)
    return ch
end
function Ψ(C, C10, D1)
    J = sqrt(det(C))
    I1 = tr(C)
    I1_bar = I1 * J^(-2 / 3)
    return C10 * (I1_bar - 3) + (1 / D1) * (J - 1)^2
end

function constitutive_driver(C, C10, D1)
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, C10, D1), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end

function make_constitutive_driver(C10, D1)
    return C -> constitutive_driver(C, C10, D1)
end

input.model_type = :threeD # or :plane_strain or :threeD
input.E , input.ν = 20.0, 0.45
E = input.E
ν = input.ν
C10 = E / (4 * (1 + ν))
D1 = 6.0 * (1.0 - 2.0 * ν) / E
input.material = make_constitutive_driver(C10, D1)

Lx, Ly, Lz = 3.17, 1.19, 0.53
nx, ny, nz = 15, 15, 6
grid = create_grid(Lx, Ly, Lz, nx, ny, nz)

input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh )
# Create CellValues and FacetValues
input.cell_values, input.facet_values = create_values()

input.ΓN = getfacetset(grid, "right")
input.facetsets = [input.ΓN]
input.traction = [.2, 0.0, 0.0]
input.tractions = Dict(1 => input.traction)

input.tol = 1e-6
input.n_load_steps = 4
input.n_iter_NR = 20
input.filename = "3D_Hyper"
input.output_dir= "/Users/aminalibakhshi/Desktop/vtu_geo/"

sol = run_fem(input);


U = sol.U_steps[end]

# Split displacements into x, y, z components
ux = U[1:3:end]
uy = U[2:3:end]
uz = U[3:3:end]

# Print max deformation
@info "Max ux = $(maximum(ux))"
@info "Max uy = $(maximum(uy))"
@info "Max uz = $(maximum(uz))"

# Print min deformation
@info "Min ux = $(minimum(ux))"
@info "Min uy = $(minimum(uy))"
@info "Min uz = $(minimum(uz))"


