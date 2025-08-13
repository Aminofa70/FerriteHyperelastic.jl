using Revise
using Ferrite
using FerriteHyperelastic
using GLMakie
using GeometryBasics
using Colors
function create_grid(Lx, Ly, Lz, nx, ny, nz)
    left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
    right = Ferrite.Vec{3}((Lx, Ly, Lz))
    grid = generate_grid(Hexahedron, (nx, ny, nz), left, right)
    addnodeset!(grid, "clamped", x -> x[2] ≈ 0.0)  # optional
    return grid
end

# Function to create cell and facet values
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

# Function to create a DOF handler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefHexahedron,1}()^3)
    Ferrite.close!(dh)
    return dh
end


Lx, Ly, Lz = 10.0, 10.0, 10.0
nx, ny, nz = 5, 5, 5
grid = create_grid(Lx, Ly, Lz, nx, ny, nz)

V, F = FerriteHyperelastic.to_geometry(grid, Ferrite.Hexahedron)

# Visualisation

fig = Figure(size=(1000, 800))
ax = Axis3(fig[1, 1],
    aspect=:data,
    xlabel="X",
    ylabel="Y",
    zlabel="Z",
    title="3D Hexahedral Mesh",
    limits=(-1, 12, -1, 12, -1, 12))


poly!(ax, GeometryBasics.Mesh(V, F),
    color=:gray,
    strokecolor=:black,
    strokewidth=3.0,
    shading=true,
    transparency=false)


############## 
outer_faces = get_boundary_faces(grid, Ferrite.Hexahedron)

ax2 = Axis3(fig[1, 2],
    aspect=:data,
    xlabel="X",
    ylabel="Y",
    zlabel="Z",
    title="Boundary condition(node)",
    limits=(-1, 12, -1, 12, -1, 12))


poly!(ax2, GeometryBasics.Mesh(V, outer_faces);
    color=(Gray(0.95), 0.3),   # light gray, 30% opacity (Gray(0.95), 0.3)
    strokecolor=:black,        # keep mesh edges visible
    strokewidth=2.0,
    shading=true,
    transparency=true          # important for alpha blending
)

## nodeset
dh = create_dofhandler(grid)
boundary = getnodeset(dh.grid, "clamped")

nodesset = to_boundary(grid, boundary, Nodes, Ferrite.Hexahedron)
# Plot clamped nodes on the surface
scatter!(ax2, nodesset,
    color=:blue,
    markersize=15.0,
    marker=:circle,
    strokecolor=:black,
    strokewidth=2,
    label="Clamped Nodes")

facets = getfacetset(grid, "right")
facessset = to_boundary(grid, facets, Faces, Ferrite.Hexahedron)

axislegend(ax2, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(fig)
