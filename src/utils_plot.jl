struct Faces end
struct Nodes end
"""
    to_geometry(grid, Ferrite.Triangle)

Convert a Ferrite triangular grid to GeometryBasics-compatible vertices and faces for plotting.
"""
function to_geometry(grid, ::Type{Ferrite.Triangle})
    ferrite_cells = grid.cells
    ferrite_nodes = grid.nodes
    F = [TriangleFace{Int64}(t.nodes...) for t in ferrite_cells]
    V = [Point{2, Float64}(Ferrite.get_node_coordinate(n)[1:2]...) for n in ferrite_nodes]
    return V, F
end

####################### end of to_geometry: Trianlge #########################
"""
    to_geometry(grid,Ferrite.Quadrilateral) 

Convert a Ferrite Quadrilateral grid to GeometryBasics-compatible vertices and faces for plotting.
"""
function to_geometry(grid, ::Type{Ferrite.Quadrilateral})
    ferrite_cells = grid.cells
    ferrite_nodes = grid.nodes
    F = [QuadFace{Int64}(t.nodes...) for t in ferrite_cells]
    V = [Point{2, Float64}(Ferrite.get_node_coordinate(n)[1:2]...) for n in ferrite_nodes]
    return V, F
end

####################### end of to_geometry: Quadrilateral #########################
"""
      to_boundary(grid, facets, Faces,Ferrite.Quadrilateral)
or 
      to_boundary(grid, facets, Faces,Ferrite.Triangle)

Convert boundary facets from a Ferrite grid to 2D point coordinates for plotting.
"""
function to_boundary(grid, facets, ::Type{Faces},::Type{T}) where {
    T<:Union{Ferrite.Quadrilateral, Ferrite.Triangle}
}
    facet_points = Point2f[]
    for facet in facets
        cell = grid.cells[facet[1]]
        facet_nodes = Ferrite.facets(cell)[facet[2]]
        for n in facet_nodes
            push!(facet_points, Point2f(Ferrite.get_node_coordinate(grid.nodes[n]).data...))
        end
    end
    return facet_points
end

## if want segement use: (not the case here)
# function to_boundary(grid, facets, ::Type{Faces})
#     segments = Vector{Vector{Point2f}}(undef, length(facets))
#     for (i, facet) in enumerate(facets)
#         cell = grid.cells[facet[1]]
#         facet_nodes = Ferrite.facets(cell)[facet[2]]
#         p1 = Point2f(Ferrite.get_node_coordinate(grid.nodes[facet_nodes[1]]).data...)
#         p2 = Point2f(Ferrite.get_node_coordinate(grid.nodes[facet_nodes[2]]).data...)
#         segments[i] = [p1, p2]
#     end
#     return segments
# end
#########  end of  to_boundary(grid, facets, ::Type{face})
"""
    to_boundary(grid, nodeset, Nodes, Ferrite.Quadrilateral )

    or

    to_boundary(grid, nodeset, Nodes, Ferrite.Triangle )

Convert a set of node indices from a Ferrite grid to 2D point coordinates for plotting.
"""
function to_boundary(grid, nodeset, ::Type{Nodes}, ::Type{T}) where {
    T<:Union{Ferrite.Quadrilateral, Ferrite.Triangle}
}
    nodesset = [Point2f(Tuple(Ferrite.get_node_coordinate(grid.nodes[i]))) for i in nodeset]

    return nodesset
end

##############################################################################
####################################### three dimensional functions

"""
     get_faces(cell)

function to extract all quadrilateral faces from hexahedral cells
"""
function get_faces(cell::Hexahedron)
    faces = [
        (cell.nodes[1], cell.nodes[2], cell.nodes[3], cell.nodes[4]),  # -z face
        (cell.nodes[5], cell.nodes[6], cell.nodes[7], cell.nodes[8]),  # +z face
        (cell.nodes[1], cell.nodes[2], cell.nodes[6], cell.nodes[5]),  # -y face
        (cell.nodes[3], cell.nodes[4], cell.nodes[8], cell.nodes[7]),  # +y face
        (cell.nodes[1], cell.nodes[4], cell.nodes[8], cell.nodes[5]),  # -x face
        (cell.nodes[2], cell.nodes[3], cell.nodes[7], cell.nodes[6])   # +x face
    ]
    return [QuadFace{Int}(f...) for f in faces]
end
##################################################
"""
     to_geometry(grid, Ferrite.Hexahedron)

Convert a Ferrite Hexahedron grid to GeometryBasics-compatible vertices and faces for plotting.
"""
function to_geometry(grid, ::Type{Ferrite.Hexahedron})
    
    # Collect all unique faces
    all_faces = QuadFace{Int}[]
    for cell in grid.cells
        append!(all_faces, get_faces(cell))
    end
    F = unique(all_faces)
    V = [Point{3,Float64}(Ferrite.get_node_coordinate(n)...) for n in grid.nodes]
    return V, F
end
##################################################
"""
    to_boundary(grid, nodeset, Nodes, Ferrite.Hexahedron) 
    or 
    to_boundary(grid, nodeset, Nodes, Ferrite.Tetrahedron) 

Convert a set of node indices from a Ferrite grid to 3D point coordinates for plotting.
"""
function to_boundary(grid, nodeset, ::Type{Nodes}, ::Type{T}) where {
    T<:Union{Ferrite.Hexahedron, Ferrite.Tetrahedron}
}
    nodesset = [Point{3, Float64}(Ferrite.get_node_coordinate(grid.nodes[i]).data...) for i in nodeset]

    return nodesset
end
##################################################
"""
       canonical_quad(t)

Canonicalize a quad (4 node indices) so that rotations/reversals map to the same key.
This ensures internal faces (shared by 2 cells) are counted correctly.
"""
function canonical_quad(t::NTuple{4,Int})
    a1 = t
    a2 = (t[2], t[3], t[4], t[1])
    a3 = (t[3], t[4], t[1], t[2])
    a4 = (t[4], t[1], t[2], t[3])
    r  = (t[4], t[3], t[2], t[1])
    r2 = (r[2], r[3], r[4], r[1])
    r3 = (r[3], r[4], r[1], r[2])
    r4 = (r[4], r[1], r[2], r[3])
    return minimum((a1,a2,a3,a4,r,r2,r3,r4))
end
##################################################
"""
    get_boundary_faces(grid , Ferrite.Hexahedron)

Function to get faces at the boundary for Ferrite.Hexahedron 
"""
function get_boundary_faces(grid,  ::Type{Ferrite.Hexahedron})
    face_counts = Dict{NTuple{4,Int},Int}()
    rep_face    = Dict{NTuple{4,Int},NTuple{4,Int}}()

    for cell in grid.cells
        for f in get_faces(cell)
            key = canonical_quad(Tuple(f))  # <- works for tuples *and* QuadFace
            face_counts[key] = get(face_counts, key, 0) + 1
            if !haskey(rep_face, key)
                rep_face[key] = Tuple(f)    # store a tuple rep
            end
        end
    end

    outer_keys  = [k for (k, c) in face_counts if c == 1]
    outer_faces = [QuadFace{Int}(rep_face[k]...) for k in outer_keys]
    return outer_faces
end
##################################################
"""
    to_boundary(grid, facets, Faces,Ferrite.Hexahedron)
    or
    to_boundary(grid, facets, Faces, Ferrite.Tetrahedron)


Convert boundary facets from a Ferrite grid to 3D point coordinates for plotting.
"""
function to_boundary(grid, facets, ::Type{Faces}, ::Type{T}) where {
    T<:Union{Ferrite.Hexahedron,Ferrite.Tetrahedron}
}
    facet_points = Point3f[]
    for facet in facets
        cell = grid.cells[facet[1]]
        facet_nodes = Ferrite.facets(cell)[facet[2]]
        for n in facet_nodes
            push!(facet_points, Point3f(Ferrite.get_node_coordinate(grid.nodes[n]).data...))
        end
    end
    facet_points = unique(facet_points)
    return facet_points
end
##################################################
"""
tet_faces_as_tuples(cell)

function to get faces for the tet4 element 
"""
function tet_faces_as_tuples(cell::Ferrite.Tetrahedron)
    n = cell.nodes
    raw = [
        (n[1], n[2], n[3]),
        (n[1], n[2], n[4]),
        (n[1], n[3], n[4]),
        (n[2], n[3], n[4]),
    ]
    return [Tuple(sort(collect(t))) for t in raw]
end
##################################################
"""
    to_geometry(grid, Ferrite.Tetrahedron)

Convert a Ferrite Hexahedron grid to GeometryBasics-compatible vertices and faces for plotting.
"""
function to_geometry(grid, ::Type{Ferrite.Tetrahedron})

    all_face_tuples = Tuple{Int,Int,Int}[]
    for cell in grid.cells
        append!(all_face_tuples, tet_faces_as_tuples(cell))
    end
    unique_face_tuples = unique(all_face_tuples)
    F = [TriangleFace{Int}(t...) for t in unique_face_tuples]
    V = [Point{3,Float64}(Ferrite.get_node_coordinate(n)...) for n in grid.nodes]
    return V, F
end
##################################################
"""
   canonical_tri(t)

# Canonical form for a triangle (handles rotations + reversal = 6 cases)
"""
function canonical_tri(t::NTuple{3,Int})
    a1 = t
    a2 = (t[2], t[3], t[1])
    a3 = (t[3], t[1], t[2])
    r  = (t[1], t[3], t[2])
    r2 = (r[2], r[3], r[1])
    r3 = (r[3], r[1], r[2])
    return minimum((a1, a2, a3, r, r2, r3))
end
##################################################
"""
    get_boundary_faces_tet(grid,  Ferrite.Tetrahedron)

Function to get faces at the boundary for Ferrite.Tetrahedron
"""
function get_boundary_faces(grid,  ::Type{Ferrite.Tetrahedron})
    face_counts = Dict{NTuple{3,Int}, Int}()
    rep_face    = Dict{NTuple{3,Int}, NTuple{3,Int}}()

    for cell in grid.cells
        for f in tet_faces_as_tuples(cell)
            key = canonical_tri(Tuple(f))          # works for tuples & TriangleFace
            face_counts[key] = get(face_counts, key, 0) + 1
            if !haskey(rep_face, key)
                rep_face[key] = Tuple(f)
            end
        end
    end

    outer_keys  = [k for (k, c) in face_counts if c == 1]
    outer_faces = [TriangleFace{Int}(rep_face[k]...) for k in outer_keys]
    return outer_faces
end
