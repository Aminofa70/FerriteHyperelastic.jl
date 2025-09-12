using Revise 
using Ferrite
using FerriteHyperelastic
using GLMakie
using GeometryBasics
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)),
        Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)),
        Ferrite.Vec{2}((0.0, Ly))
    ]
    Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
end

grid = create_grid(10, 10,30, 30)

V, F = FerriteHyperelastic.to_geometry(grid, Ferrite.Quadrilateral)

# Visualisation
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(),
    xlabel="X", ylabel="Y",
    title="Mesh with Boundary Conditions",
    limits=(-1, 12, -1, 12))


poly!(ax, GeometryBasics.Mesh(V, F),
    color=:gray,
    strokecolor=:black,
    strokewidth=1,
    shading=false)


# add boundary condition with nodes
addnodeset!(grid, "clamped", x -> x[1] ≈ 0.0)

clamped_nodes = getnodeset(grid, "clamped")

nodesset = FerriteHyperelastic.to_boundary(grid, clamped_nodes, Nodes, Ferrite.Quadrilateral)
scatter!(ax, nodesset,
    color=:red,
    markersize=10,
    marker=:rect,
    strokecolor=:black,
    strokewidth=2,
    label="Clamped Nodes",)
axislegend(ax, position=:lt, backgroundcolor=(:white, 0.7), framecolor=:gray)

display(fig)
