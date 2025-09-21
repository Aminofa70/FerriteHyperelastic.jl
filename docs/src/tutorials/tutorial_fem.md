## Guid to use Finite Element (fem)

### Create Input Data Structure
Create a dynamic structure for putting FEM parameters into it.
```
input = InputStruct()
```
!!! tipo
    Define geometry, mesh, and FEM parameters for interpolation and numerical integration, and define degree-of-freedom using [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/).

    More details for these parts can be found in [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/).

### Mesh generation
Generate mesh 
```
"""
Lx : Lenght in x-direction
Ly : Lenght in y-direction
Ferrite.Quadrilateral: Define the type of elements.
nx : Number of element in x-dir 
ny : Number of element in y-dir
"""
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

```

### Define FEM values 
Define FEM values (interpolation and numerical integration)
```
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral,order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(2)
    return Ferrite.CellValues(qr, ip), Ferrite.FacetValues(qr_face, ip)
end

```
!!! note
    Very important to know is that the type here should be the same as the element type.

### Define DOF
Define dof (degree of freedom)
```
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral,1}()^2)
    Ferrite.close!(dh)
    return dh
end

```
### Dirichlet Boundary Condition 
Define Dirichlet Boundary Condition

```
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "support_1"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.close!(ch)
    return ch
end

```
!!! hint
    Optional step: if users need to plot a force-displacement curve, getting the degree of freedom from reaction force and displacement is required. For the displacement, the Dirichlet boundary condition is used. For the force, another boundary condition is defined. If this plot is not the aim, it is not necessary to define the following code. 

```
function create_bc_force(dh)
    dbc = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(dbc, Ferrite.Dirichlet(:u, getfacetset(grid, "pressure"), (x, t) -> 0*x))
    Ferrite.close!(dbc)
    return dbc
end

```

!!! important
    It seems that Ferrite.jl does not have a direct way for finding the degree of freedom. For this reason, the above code is defined to find the degree of freedom of the reaction force, as will be seen in the following.

To now, the required functions for FEM were defined.

These parts are from the package [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/). 

***
***
### Define hyperelastic strain energy function (here neo-Hookean)
The strain energy function is defined
```
function Ψ(C, C10, D1)
    J = sqrt(det(C))
    I1 = tr(C)
    I1_bar = I1 * J^(-2 / 3)
    return C10 * (I1_bar - 3) + (1 / D1) * (J - 1)^2
end

```
Second Piola–Kirchhoff stress tensor
```
function constitutive_driver(C, C10, D1)
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, C10, D1), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end

```
And for calling the strain energy using C, the first invariant of the Cauchy-Green deformation tensor in the FEM solver

```
function make_constitutive_driver(C10, D1)
    return C -> constitutive_driver(C, C10, D1)
end

```
!!! hint 
    In the following, the FEM parameters are put into the dynamic structure defined above. 

### Problem and Load Type
Problem and load type are defined

!!! note "Problem Types"
    *  :plane_strain
    *  :plane_stress
    *  :threeD

!!! note "Load types"
    * :traction
    * :displacement

This example is ```:plane_strain```  and  ``` :traction ```.

```
input.model_type = :plane_strain   
input.load_type = :traction

```

#### Define material parameters and Stress
```
# E: Young's Modulus
# ν: Poisson's ratio

input.E , input.ν = 3.35, 0.45
E = input.E
ν = input.ν

C10 = E / (4 * (1 + ν))
D1 = 6.0 * (1.0 - 2.0 * ν) / E
input.material = make_constitutive_driver(C10, D1)

```

#### Generate Geometry and Mesh
Call the mehs function and put into in the dynamic structure.
```
# Define parameters for the plate and mesh
Lx, Ly = 3.17, 1.73  # Plate dimensions
nx, ny = 10, 10   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid

```
#### FEM Values, DOF, Dirichlet BC
Get FEM values, DOF, Dirichlet BC and put them into the dynamic structure.
```
input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh )
# Create CellValues and FacetValues
input.cell_values, input.facet_values = create_values()

```
#### Traction Force
```
input.ΓN = getfacetset(grid, "pressure")
input.facetsets = [input.ΓN]
input.traction = [2.2, 0.0]
input.tractions = Dict(1 => input.traction)

```
#### Dofs Reaction Force & Displacement
To plot force-displacement, Dofs Reaction Force & Displacement are put in the dynamic structure
```
dof_F_x = input.ch.prescribed_dofs[1:2:end]
input.dof_F = dof_F_x;


dbc= create_bc_force(input.dh)
dof_U_x = dbc.prescribed_dofs[1:2:end] 
input.dof_U = dof_U_x

```
!!! note
    If the force-displacement plot is not required, then

    ```
    input.dof_F = []
    input.dof_U = []
    ```
#### Solver Parameters
The settings for the solver is defined as
```
input.tol = 1e-6
maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver(500,1.0,1e-3,1e-15,0.8,1000)

input.maxIterPerInc = maxIterPerInc
input.totalTime = totalTime
input.initInc = initInc
input.minInc = minInc
input.maxInc = maxInc
input.totalInc = totalInc

```
!!! hint
    The default value of the solver can also be used as

    ```
    maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver()
    ```
#### Save VTU file
The displacement of the last step is saved in the VTU file. Define a name and directory for the file. However, the built-in plot can also be implemented. 
```
input.filename = "2D_Hyper"
input.output_dir= "/Users/aminalibakhshi/Desktop/vtu_geo/"
```

#### Run FEM solver
```
sol  = run_fem(input)
```
!!! note 
    The solver returns a structure with the following fields:

    ```
    # fieldnames(typeof(sol)) ->(:U_steps, :U_effect, :F_effect) 

    ```
#### Find MAX Displacement

```
U = sol.U_steps[end]

# Extract final displacement and evaluate at grid nodes
U = sol.U_steps[end]
u_nodes = vec(evaluate_at_grid_nodes(input.dh, U, :u))
ux, uy = getindex.(u_nodes, 1), getindex.(u_nodes, 2)
@info "Max |ux| = $(maximum(abs, ux)), Max |uy| = $(maximum(abs, uy))"

```
 

#### Plot Displacement Force

```
GLMakie.closeall()
fig = Figure(size=(800, 600), fontsize=26)
ax = Axis(fig[1, 1], xlabel="Displacement", ylabel="Force", title="force-displacement", xgridvisible = false, ygridvisible = false)
lines!((sol.U_effect), abs.(sol.F_effect), color = :black)
scatter!((sol.U_effect), abs.(sol.F_effect), marker = :circle , color = :red)
display(fig)

```