using Revise
using FerriteHyperelastic
using Ferrite
using GLMakie
using GeometryBasics
using ColorSchemes
const Vec = Ferrite.Vec

##################################################################e
input = InputStruct()
##################################################################
function create_cook_grid(nx, ny)
    corners = [
        Vec{2}((0.0, 0.0)),
        Vec{2}((48.0, 44.0)),
        Vec{2}((48.0, 60.0)),
        Vec{2}((0.0, 44.0)),
    ]
    grid = generate_grid(Quadrilateral, (nx, ny), corners)
    # facesets for boundary conditions
    addfacetset!(grid, "clamped", x -> norm(x[1]) ≈ 0.0)
    addfacetset!(grid, "traction", x -> norm(x[1]) ≈ 48.0)
    return grid
end;
##################################################################
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral,order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(1)
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
    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getfacetset(dh.grid, "clamped"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.close!(ch)
    return ch
end
##################################################################
function Ψ(C, C10, D1)
    J = sqrt(det(C))
    if det(C) < 0
        error("determinant is negative")
    end
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
input.model_type = :plane_stress  # or :plane_strain or ::threeD
input.load_type = :traction

input.E, input.ν = 4.35, 0.45
E = input.E
ν = input.ν
C10 = E / (4 * (1 + ν))
D1 = 6.0 * (1.0 - 2.0 * ν) / E
input.material = make_constitutive_driver(C10, D1)
##################################################################

nx, ny = 10, 10 # Number of elements along x and y
grid = create_cook_grid(nx, ny)

input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh)
# Create CellValues and FacetValues
input.cell_values, input.facet_values = create_values()

input.ΓN = getfacetset(grid, "traction")
input.facetsets = [input.ΓN]
input.traction = [0.0, 1.17]
input.tractions = Dict(1 => input.traction)
##################################################################
##################################

input.dof_F = []
input.dof_U = []
##################################################################
input.tol = 1e-6

## default
maxIterPerInc, totalTime, initInc, minInc, maxInc, totalInc = initialize_solver()

# change like the following if you need
# maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver(500,1.0,1e-3,1e-15,0.2,1000)

input.maxIterPerInc = maxIterPerInc
input.totalTime = totalTime
input.initInc = initInc
input.minInc = minInc
input.maxInc = maxInc
input.totalInc = totalInc
##################################################################

input.filename = "2D_Hyper"
input.output_dir = "/Users/aminalibakhshi/Desktop/vtu_geo/"
##################################################################
################  solution 
sol = run_fem(input)

V, F = FerriteHyperelastic.to_geometry(grid, Ferrite.Quadrilateral)


function plot_disp()
    GLMakie.closeall()
    fig = Figure(size = (1000, 600))
    
    ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "X", ylabel = "Y")
    
    poly!(ax, GeometryBasics.Mesh(V, F), 
          color = (:gray, 0.10), strokecolor = :black, strokewidth = 1, shading = false)
    
    incRange = length(sol.U_steps) 
    
    UT        = fill(V, incRange)                      
    UT_mag    = fill(zeros(length(V)), incRange)
    ut_mag_max = zeros(incRange)                       
    
    @inbounds for i in 1:incRange
        U = sol.U_steps[i]
        u_nodes = vec(evaluate_at_grid_nodes(input.dh, U, :u))  
        ux = getindex.(u_nodes, 1)
        uy = getindex.(u_nodes, 2)
        DD_disp = Vector{Vector{Float64}}()
        for j in eachindex(ux)
            push!(DD_disp, [ux[j], uy[j]])
        end
        UT[i] = [Point{2,Float64}(u) for u in DD_disp]    
        UT_mag[i] = norm.(UT[i])                         
        ut_mag_max[i] = maximum(UT_mag[i])              
    end
    
    
    scale = 0.0  
    
    step_index = Observable(1)
    Vdef_obs = Observable(V .+ scale .* UT[1])
    color_obs = Observable(UT_mag[1])
    colorrange_obs = Observable((0.0, ut_mag_max[1]))
    
    
    mobj = poly!(ax, 
     lift(Vdef_obs) do verts
            GeometryBasics.Mesh(verts, F)
        end,
        color = color_obs,
        strokewidth = 2,
        transparency = false,
        colormap = Reverse(:Spectral),
        colorrange = colorrange_obs
    )
    
    Colorbar(fig[1, 2], mobj, label = "Displacement magnitude [mm]")
    
    slider = Slider(fig[2, 1], range = 1:incRange, startvalue = 1, width = 800, linewidth=30)
    step_label = Label(fig[2, 2], "Step: 1")
    
    on(slider.value) do i
        step_index[] = i
        Vdef_obs[] = V .+ scale .* UT[i]
        color_obs[] = UT_mag[i]
        colorrange_obs[] = (0.0, ut_mag_max[i])  
        step_label.text = "Step: $i"
    end
    display(fig) 
end

plot_disp()

