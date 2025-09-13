# This is the instruction to use the package
## Guid to use Finite Element(fem)
First we need to activate a dynamic structure to have all finite element parameters in it. We do it using 
```
input = InputStruct()
```

Here ```input``` is the a structure that will contains all fem parameters. 

Now using [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/) we define geometry, mesh, fem parameters for interpolation and numerical integration and define degree-of-freedom. More details for this parts can be found in [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/). 

As an example, we do it for a two dimensional problem.
Our element here is quad4 and then we define mesh and geometry as
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
It is noted that 
```
addnodeset!(grid, "support_1", x -> x[1] ≈ 0.0)
addfacetset!(grid, "pressure", x -> x[1] ≈ Lx)

```
are used to define our boundary conditions.

If you wan to chose other element types like tri3 and so on please refer to [generate_grid](https://github.com/Ferrite-FEM/Ferrite.jl/blob/f1d1d0deef7bdaf019bd63ce9e8d959b6ebc8c4d/src/Grid/grid_generators.jl#L1-L7) .

Next step is defining the fem values (interpolation and numerical integration)

```
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral,order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(2)
    return Ferrite.CellValues(qr, ip), Ferrite.FacetValues(qr_face, ip)
end
```
Very important to know is that the type here should be the same as the element type. For example here we have
``` Ferrite.Quadrilateral ``` and ``` Ferrite.RefQuadrilateral ``` .

Then defining the dof (here is 2D, ```ux```,```uy```)

```
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral,1}()^2)
    Ferrite.close!(dh)
    return dh
end

```
Now, we define the Dirichlet boundary condition 

```
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, "support_1"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.close!(ch)
    return ch
end

```

If we want to plot force-displacement curve, we need to define another boundary condition to get dof of the displacment for the plot.
Displacement where the force is applied. 
```
function create_bc_force(dh)
    dbc = Ferrite.ConstraintHandler(dh)
    Ferrite.add!(dbc, Ferrite.Dirichlet(:u, getfacetset(grid, "pressure"), (x, t) -> 0*x))
    Ferrite.close!(dbc)
    return dbc
end

```
To now we defined the required functions for fem. These parts are from the package [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/). More details can be seen in this package. 

***

The hyperelastic strain enenrgy function and stress are now defined.
We use the neo-Hookean function, but other strain energy functions can also be employed

```
function Ψ(C, C10, D1)
    J = sqrt(det(C))
    I1 = tr(C)
    I1_bar = I1 * J^(-2 / 3)
    return C10 * (I1_bar - 3) + (1 / D1) * (J - 1)^2
end

```

The second Piola kirchhoff stress and its differentiation are then define 

```
function constitutive_driver(C, C10, D1)
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, C10, D1), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end

```
Now, we define a driver for the strain energy in order to inlcude it in the input structure

```
function make_constitutive_driver(C10, D1)
    return C -> constitutive_driver(C, C10, D1)
end
```

We should tell the solve our problem is 2D or 3D. Also we should define we have applied traction or displacement. Then we should script 

```
input.model_type = :plane_strain   # or :plane_strain; :plane_stress; :threeD
input.load_type = :traction

```

Now we should assign parameters for the materials ans also call the fem function to include them in the input structure.

```

input.E , input.ν = 3.35, 0.45
E = input.E
ν = input.ν
C10 = E / (4 * (1 + ν))
D1 = 6.0 * (1.0 - 2.0 * ν) / E
input.material = make_constitutive_driver(C10, D1)


Lx, Ly = 3.17, 1.73  # Plate dimensions
nx, ny = 10, 10   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid

input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh )
# Create CellValues and FacetValues
input.cell_values, input.facet_values = create_values()

```
Apply the traction (it can be multiple tractions)

```
input.ΓN = getfacetset(grid, "pressure")
input.facetsets = [input.ΓN]
input.traction = [2.2, 0.0]
input.tractions = Dict(1 => input.traction)

```
If the aim is also plot force-displacement, we shoul also find the dof of reaction force and displacement (if not leave them empty)

```

dof_F_x = input.ch.prescribed_dofs[1:2:end]
input.dof_F = dof_F_x;

dbc= create_bc_force(input.dh)
dof_U_x = dbc.prescribed_dofs[1:2:end] 
input.dof_U = dof_U_x
# input.dof_F = []
# input.dof_U = []

```
Define the tolernce for solver and also parameters for time integration. Because the problem is nonlinear we need to solve it incrementally (time step here), then
```
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

```
The solver saves the displacement of the last time step in vtu file for analysis of the results. So solver needs a name and dir to save it 
```
input.filename = "2D_Hyper"
input.output_dir= "/Users/aminalibakhshi/Desktop/vtu_geo/"

```
Howover, the solver saves the results for each time steps and returns all corresponding displacement vector and use can plot it using Plot.jl or GLMakie.jl

Now, calling the solve

```
sol  = run_fem(input)

```

Solver is a structure that returns ``` U = sol.U_steps[end] ``` the displacement for each step and ``` sol.F_effect ``` the sum of reaction force and ```sol.U_effect ``` average displacement in applied force (or displacement).

The force-displacement can be plotted using GLMakie.jl as

```
GLMakie.closeall()
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], xlabel="Displacement", ylabel="Force", title="force-displacement", xgridvisible = false, ygridvisible = false)



lines!((sol.U_effect), abs.(sol.F_effect), color = :black)
scatter!((sol.U_effect), abs.(sol.F_effect), marker = :circle , color = :red)
display(fig)

```







