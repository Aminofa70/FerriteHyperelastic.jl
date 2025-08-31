using Revise
using FerriteHyperelastic
using Ferrite
##################################################################
input = InputStruct()
##################################################################
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "clamped", x -> x[1] ≈ 0.0)
    addfacetset!(grid, "shear", x -> x[1] ≈ Lx)
    addfacetset!(grid, "bending", x -> x[2] ≈ Ly)
    return grid
end
##################################################################
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral,order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(2)
    return Ferrite.CellValues(qr, ip), Ferrite.FacetValues(qr_face, ip)
end
##################################################################
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral,1}()^2)
    Ferrite.close!(dh)
    return dh
end
##################################################################
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "clamped"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.close!(ch)
    return ch
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
##################################################################
input.model_type = :plane_strain   # or :plane_strain or :threeD
input.load_type = :traction
##################################################################
input.E , input.ν = 4.35, 0.49
E = input.E
ν = input.ν
C10 = E / (4 * (1 + ν))
D1 = 6.0 * (1.0 - 2.0 * ν) / E
input.material = make_constitutive_driver(C10, D1)
##################################################################

# Define parameters for the plate and mesh
Lx, Ly = 3.0, 1.0  # Plate dimensions
nx, ny = 30, 30   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid

input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh )
# Create CellValues and FacetValues
input.cell_values, input.facet_values = create_values()

input.ΓN = getfacetset(grid, "shear")
input.ΓD = getfacetset(grid, "bending")
input.facetsets = [input.ΓN , input.ΓD ]
input.bending = [0.0, -0.05]
input.shear= [0.0, -0.1]
input.tractions = Dict(
    1 => input.shear,
    2 => input.bending)
##################################################################

input.tol = 1e-6

## default
#maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver()

# change like the following if you need
maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver(500,1.0,1e-3,1e-15,0.2,1000)
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
sol = run_fem(input)

U = sol.U_steps[end]
# Split displacements into x and y components
ux = U[1:2:end]
uy = U[2:2:end]

# Print max deformation if desired
@info "Max ux = $(maximum(ux))"
@info "Max uy = $(maximum(uy))"

# Print min deformation if desired
@info "Min ux = $(minimum(ux))"
@info "Min uy = $(minimum(uy))"



