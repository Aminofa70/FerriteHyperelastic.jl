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
    grid = Ferrite.generate_grid(Ferrite.Triangle, (nx, ny), corners)
    return grid
end


grid = create_grid(10, 10, 5, 5)


V, F = FerriteHyperelastic.to_geometry(grid, Ferrite.Triangle)

# Visualisation
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(), 
          xlabel="X", ylabel="Y", 
          title="Mesh with Boundary Conditions",
          limits=(-1, 12, -1, 12))


poly!(ax, GeometryBasics.Mesh(V, F), 
      color=:lightblue, 
      strokecolor=:black,
      strokewidth=1,
      shading=false)


left_facets = getfacetset(grid, "top")
facet_points = FerriteHyperelastic.to_boundary(grid, left_facets, Faces, Ferrite.Triangle)

scatter!(ax, facet_points;
    color=:black,
    markersize=15,
    marker=:circle,
    strokecolor=:black,
    strokewidth=2,
    label="top Facet Nodes"
)
axislegend(ax, position=:lt, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(fig)