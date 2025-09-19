using Revise
using FerriteHyperelastic
using Ferrite
using GLMakie
using GeometryBasics 
##################################################################
# create the structure for saving fem input 
input = InputStruct()
##################################################################
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "support_1", x -> x[1] ≈ 0.0)
    addfacetset!(grid, "pressure", x -> x[1] ≈ Lx)
    return grid
end
##################################################################
# create fem values: interpolation and quadrature point
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral,order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(2)
    return Ferrite.CellValues(qr, ip), Ferrite.FacetValues(qr_face, ip)
end
##################################################################
# create dof (degree of freedom)
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral,1}()^2)
    Ferrite.close!(dh)
    return dh
end
##################################################################
# creat boundary condition Dirichlet (displacement)
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "support_1"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.close!(ch)
    return ch
end

####################################
# creat a seperate boundary for load to get  dof
# creat boundary condition Dirichlet (displacement)
function create_bc_force(dh)
    dbc = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(dbc, Ferrite.Dirichlet(:u, getfacetset(grid, "pressure"), (x, t) -> 0*x))
    Ferrite.close!(dbc)
    return dbc
end
##################################################################
function Ψ(C, C10, D1)
    J = sqrt(det(C))
    I1 = tr(C)
    I1_bar = I1 * J^(-2 / 3)
    return C10 * (I1_bar - 3) + (1 / D1) * (J - 1)^2
end
##################################################################
function constitutive_driver(C, C10, D1)
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, C10, D1), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end
##################################################################
function make_constitutive_driver(C10, D1)
    return C -> constitutive_driver(C, C10, D1)
end

input.model_type = :plane_strain   # or :plane_strain; :plane_stress; :threeD
input.load_type = :traction
##################################################################

input.E , input.ν = 3.35, 0.45
E = input.E
ν = input.ν

C10 = E / (4 * (1 + ν))
D1 = 6.0 * (1.0 - 2.0 * ν) / E
input.material = make_constitutive_driver(C10, D1)
##################################################################

# Define parameters for the plate and mesh
Lx, Ly = 3.17, 1.73  # Plate dimensions
nx, ny = 10, 10   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid

input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh )
# Create CellValues and FacetValues
input.cell_values, input.facet_values = create_values()

input.ΓN = getfacetset(grid, "pressure")
input.facetsets = [input.ΓN]
input.traction = [2.2, 0.0]
input.tractions = Dict(1 => input.traction)

##################################
dof_F_x = input.ch.prescribed_dofs[1:2:end]
input.dof_F = dof_F_x;


dbc= create_bc_force(input.dh)
dof_U_x = dbc.prescribed_dofs[1:2:end] 
input.dof_U = dof_U_x

# input.dof_F = []
# input.dof_U = []
##################################################################
input.tol = 1e-6

## default
#maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver()

# change like the following if you need
maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver(500,1.0,1e-3,1e-15,0.8,1000)


input.maxIterPerInc = maxIterPerInc
input.totalTime = totalTime
input.initInc = initInc
input.minInc = minInc
input.maxInc = maxInc
input.totalInc = totalInc
##################################################################

input.filename = "2D_Hyper"
input.output_dir= "/Users/aminalibakhshi/Desktop/vtu_geo/"
##################################################################
################  solution 

sol  = run_fem(input)

# fieldnames(typeof(sol)) ->(:U_steps, :U_effect, :F_effect) 
U = sol.U_steps[end]

# Extract final displacement and evaluate at grid nodes
U = sol.U_steps[end]
u_nodes = vec(evaluate_at_grid_nodes(input.dh, U, :u))
ux, uy = getindex.(u_nodes, 1), getindex.(u_nodes, 2)
@info "Max |ux| = $(maximum(abs, ux)), Max |uy| = $(maximum(abs, uy))"



GLMakie.closeall()
fig = Figure(size=(800, 600), fontsize=26)
ax = Axis(fig[1, 1], xlabel="Displacement", ylabel="Force", title="force-displacement", xgridvisible = false, ygridvisible = false)



lines!((sol.U_effect), abs.(sol.F_effect), color = :black)
scatter!((sol.U_effect), abs.(sol.F_effect), marker = :circle , color = :red)
display(fig)








