using Revise, Ferrite, FerriteHyperelastic, GLMakie, GeometryBasics

# Create grid function
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

# Generate grid
grid = create_grid(10, 10, 5, 5)

# Convert to plottable format
V, F = FerriteHyperelastic.to_geometry(grid, Ferrite.Triangle)

# Create figure
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(), 
          xlabel="X", ylabel="Y", 
          title="Mesh with Boundary Conditions",
          limits=(-1, 12, -1, 12))

# Plot the mesh
poly!(ax, GeometryBasics.Mesh(V, F), 
      color=:lightblue, 
      strokecolor=:black,
      strokewidth=1,
      shading=false)

# Plot left boundary
left_facets = getfacetset(grid, "left")
segments = FerriteHyperelastic.to_boundary(grid, left_facets, Faces)

# # plot boundary as a segment
# for seg in segments
#     lines!(ax, seg, color=:red, linewidth=5)
# end

### nodes the face boundary can also be shown 

# First plot all segments without labels
for seg in segments
    scatter!(ax, seg, color=:red, markersize=20, strokecolor=:red, strokewidth=2)
end

# Then plot just one point with the label
if !isempty(segments)
    scatter!(ax, [segments[1][1]],  # Just use first point
             color=:red, markersize=20, strokecolor=:red, strokewidth=2,
             label="BC")
end

# Add legend
axislegend(ax, position=:lt, backgroundcolor=(:white, 0.7), framecolor=:gray)

display(fig)