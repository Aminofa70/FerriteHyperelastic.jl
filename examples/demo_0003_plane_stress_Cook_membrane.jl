using Revise
using FerriteHyperelastic
using Ferrite
input = InputStruct()
function create_cook_grid(nx, ny)
    corners = [
        Vec{2}((0.0, 0.0)),
        Vec{2}((48.0, 44.0)),
        Vec{2}((48.0, 60.0)),
        Vec{2}((0.0, 44.0)),
    ]
    grid = generate_grid(Quadrilateral, (nx, ny), corners)
    # facesets for boundary conditions
    addfacetset!(grid, "clamped", x -> norm(x[1]) ≈ 0.0)
    addfacetset!(grid, "traction", x -> norm(x[1]) ≈ 48.0)
    return grid
end;
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral,order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(2)
    return Ferrite.CellValues(qr, ip), Ferrite.FacetValues(qr_face, ip)
end
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral,1}()^2)
    Ferrite.close!(dh)
    return dh
end
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getfacetset(dh.grid, "clamped"), (x, t) -> [0.0, 0.0], [1, 2]))
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

input.model_type = :plane_strain   # or :plane_strain or :3d

input.E , input.ν = 4.35, 0.45
E = input.E
ν = input.ν
C10 = E / (4 * (1 + ν))
D1 = 6.0 * (1.0 - 2.0 * ν) / E
input.material = make_constitutive_driver(C10, D1)

nx, ny = 30, 30   # Number of elements along x and y
grid = create_cook_grid(nx, ny)

input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh )
# Create CellValues and FacetValues
input.cell_values, input.facet_values = create_values()

input.ΓN = getfacetset(grid, "traction")
input.facetsets = [input.ΓN]
input.traction = [0.0, .57]
input.tractions = Dict(1 => input.traction)


input.tol = 1e-6
input.n_load_steps = 10
input.n_iter_NR = 50
input.filename = "2D_Hyper"
input.output_dir= "/Users/aminalibakhshi/Desktop/vtu_geo/"

sol = run_fem(input);

U = sol.u
# Split displacements into x and y components
ux = U[1:2:end]
uy = U[2:2:end]

# Print max deformation if desired
@info "Max ux = $(maximum(abs.(ux)))"
@info "Max uy = $(maximum(abs.(uy)))"


coords = input.grid.nodes

target = Vec{2}((48.0, 60.0))
tol = 1e-6

node_index = findfirst(i -> norm(coords[i].x - target) < tol, 1:length(coords))

if node_index !== nothing
    ux_node = U[2*node_index - 1]
    uy_node = U[2*node_index]
    println("Displacement at node (48.0, 60.0): ux = $ux_node, uy = $uy_node")
else
    println("Node at (48.0, 60.0) not found!")
end
