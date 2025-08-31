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
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "support_1"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getfacetset(dh.grid, "disp"), (x, t) -> [t], [1]))
    Ferrite.close!(ch)
    Ferrite.update!(ch, 0.0)
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

input.model_type = :plane_strain   # or :plane_strain or :threeD
input.load_type = :displacement

input.E , input.ν = 4.35, 0.45
E = input.E
ν = input.ν
C10 = E / (4 * (1 + ν))
D1 = 6.0 * (1.0 - 2.0 * ν) / E
input.material = make_constitutive_driver(C10, D1)

# Define parameters for the plate and mesh
Lx, Ly = 3.17, 1.73  # Plate dimensions
nx, ny = 4, 4 # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid

input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh)
# Create CellValues and FacetValues
input.cell_values, input.facet_values = create_values()

input.displacement = 0.5
input.tol = 1e-6
# input.n_load_steps = 10
input.n_iter_NR = 500
#input.LINE_SEARCH_STEPS = 10
# input.total_time = 1.0
# input.Δt = 0.05
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
sol = run_plane_strain_disp(input)

#sol = run_fem(input)

# Final displacement vector
U = sol.U_steps[end]

# Displacement at all nodes as a Vector{Vec{2,Float64}}
u_nodes = vec(evaluate_at_grid_nodes(input.dh, U, :u))

# Split into scalar components
ux = getindex.(u_nodes, 1)
uy = getindex.(u_nodes, 2)
@info "Max |ux| = $(maximum(abs.(ux)))"
@info "Max |uy| = $(maximum(abs.(uy)))"