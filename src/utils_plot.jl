struct Faces end
struct Nodes end
"""
    to_geometry(grid, ::Type{Ferrite.Triangle}) -> (vertices, faces)

Convert a Ferrite triangular grid to GeometryBasics-compatible vertices and faces for plotting.

# Arguments
- `grid`: A Ferrite grid containing triangular elements
- `::Type{Ferrite.Triangle}`: Type annotation ensuring only triangular grids are processed

# Returns
- `vertices`: Vector of 3D points as `Point{3,Float64}`
- `faces`: Vector of triangular faces as `TriangleFace{Int64}`

"""
function to_geometry(grid, ::Type{Ferrite.Triangle})
    ferrite_cells = grid.cells
    ferrite_nodes = grid.nodes
    F = [TriangleFace{Int64}(t.nodes...) for t in ferrite_cells]
    V = [Point{3,Float64}(Ferrite.get_node_coordinate(n)..., 0.0) for n in ferrite_nodes]
    return V, F
end

####################### end of to_geometry: Trianlge #########################
"""
    to_geometry(grid, ::Type{Ferrite.Quadrilateral}) -> (vertices, faces)

Convert a Ferrite Quadrilateral grid to GeometryBasics-compatible vertices and faces for plotting.

# Arguments
- `grid`: A Ferrite grid containing triangular elements
- `::Type{Ferrite.Quadrilateral}`: Type annotation ensuring only Quadrilateral grids are processed

# Returns
- `vertices`: Vector of 3D points as `Point{3,Float64}`
- `faces`: Vector of triangular faces as `TriangleFace{Int64}`

"""
function to_geometry(grid, ::Type{Ferrite.Quadrilateral})
    ferrite_cells = grid.cells
    ferrite_nodes = grid.nodes
    F = [QuadFace{Int64}(t.nodes...) for t in ferrite_cells]
    V = [Point{3, Float64}(Ferrite.get_node_coordinate(n)..., 0.0) for n in ferrite_nodes]
    return V, F
end
####################### end of to_geometry: Quadrilateral #########################
"""
to_boundary(grid, facets, ::Type{faces}) -> Vector{Vector{Point2f}}

Convert boundary facets from a Ferrite grid to line segments for plotting.

# Arguments
- `grid`: A Ferrite grid containing the boundary facets
- `facets`: Collection of FacetIndex tuples (cell_index, local_facet_index)

# Returns
- Vector of line segments, where each segment is a pair of Point2f coordinates
or plot for points

"""
function to_boundary(grid, facets, ::Type{Faces})
    segments = Vector{Vector{Point2f}}(undef, length(facets))
    for (i, facet) in enumerate(facets)
        cell = grid.cells[facet[1]]
        facet_nodes = Ferrite.facets(cell)[facet[2]]
        p1 = Point2f(Ferrite.get_node_coordinate(grid.nodes[facet_nodes[1]]).data...)
        p2 = Point2f(Ferrite.get_node_coordinate(grid.nodes[facet_nodes[2]]).data...)
        segments[i] = [p1, p2]
    end
    return segments
end

#########  end of  to_boundary(grid, facets, ::Type{face})
"""
    to_boundary(grid, nodeset::Vector{Int}, ::Type{nodes}) -> Vector{Point2f}

Convert a set of node indices from a Ferrite grid to 2D point coordinates for plotting.

# Arguments
- `grid`: A Ferrite grid containing the nodes
- `nodeset`: Vector of node indices to convert
- `::Type{nodes}`: Type annotation specifying node conversion mode

# Returns
- Vector of `Point2f` coordinates (Float32 for efficient plotting)

# Notes
- The third argument `::Type{nodes}` is a dispatch marker that distinguishes this
  from the facet-based boundary conversion
- Returns points in Float32 format (Point2f) which is optimal for GLMakie visualization

# Example
```julia
clamped_nodes = getnodeset(grid, "clamped")
points = to_boundary(grid, clamped_nodes, nodes)
"""
function to_boundary(grid, nodeset, ::Type{Nodes} )
    return [Point2f(Ferrite.get_node_coordinate(grid.nodes[i]).data...) for i in nodeset]
end

