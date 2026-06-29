using Mammo
using Mammo.Comodo
using Mammo.Comodo.GLMakie
using Mammo.Comodo.GeometryBasics
using Mammo.Comodo.Rotations
using Random
using Statistics
using Printf
using ComodoFerrite
using Ferrite
using SparseArrays
using BlockArrays
using LinearAlgebra
# ================================================================================
# Visualization settings
# ================================================================================

GLMakie.closeall()

strokewidth = 0.5
cmap = Makie.Categorical(:Spectral)
depth_shift = Float32(0.0)

# ================================================================================
# Geometry parameters
# ================================================================================

n = 3
r = 40.0
r1 = r / 2.5
r2 = r / 7.0
w = (r1 - r2) / 20.0
h = r2 / 1.5
sy = 1.0
gravityShiftPercentage = 0.5


# ================================================================================
# Load parameters
# ================================================================================

time_steps = 25
step_size = 1.0 / time_steps
numSteps = time_steps
load_time = 1.0

tissueDensity = 1.0e-9
gravityConstant = 9.81 * 1e3

gravityVector = Tensors.Vec{3,Float64}((0.0, 0.0, gravityConstant))

# ================================================================================
# Ogden parameters
# ================================================================================

bulkModulusFactor = 100.0

c1 = 2.0e-3 / 3.0
m1 = 8.0
c2 = c1
m2 = -m1
Kbulk = bulkModulusFactor * c1

# ================================================================================
# Geometry and mesh
# ================================================================================

F1, V1, C1, C_skin = breast_surface(n, r, r1, r2, w, h, sy, gravityShiftPercentage,)
pointSpacing = mean(edgelengths(F1, V1))

v_region = mean(V1)

vol1 = 20.0 * pointSpacing^3 / (6.0 * sqrt(2.0))

stringOpt = "paAqY"

E, V, CE, Fb, Cb = tetgenmesh(F1, V1; facetmarkerlist=C1, V_regions=[v_region], region_vol=vol1, stringOpt=stringOpt,)
grid = ComodoToFerrite(E, V)

indices_chest_nodes = unique(reduce(vcat, Fb[Cb .== 2]))

addfacetset!(grid, "chest", boundary_facets(grid, indices_chest_nodes))

F = element2faces(E)

println("Number of nodes = ", length(V))
println("Number of elements = ", length(E))
println("Number of chest nodes = ", length(indices_chest_nodes))

# ================================================================================
# Boundary condition
# ================================================================================

function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "chest"), (x, t) -> Tensors.Vec{3,Float64}((0.0, 0.0, 0.0)), [1, 2, 3]))
    Ferrite.close!(ch)
    Ferrite.update!(ch, 0.0)
    return ch
end

# ================================================================================
# Ogden material
# ================================================================================

struct Ogden
    c1::Float64
    c2::Float64
    m1::Float64
    m2::Float64
    K::Float64
end

function tensor2_from_matrix(A::AbstractMatrix{<:Real})
    return Tensor{2,3}(Float64[
        A[1, 1] A[1, 2] A[1, 3]
        A[2, 1] A[2, 2] A[2, 3]
        A[3, 1] A[3, 2] A[3, 3]
    ])
end

function tensor4_from_array(A::Array{Float64,4})
    return Tensor{4,3}(A)
end

function basis_tensor_2(a::Int, b::Int)
    A = zeros(3, 3)
    A[a, b] = 1.0
    return tensor2_from_matrix(A)
end

function finite_tensor(x)
    return all(isfinite, collect(x))
end

function has_repeated_eigs(F; tol=1.0e-8)
    C = tdot(F)
    eigC = eigvals(Symmetric(Matrix(C)))

    if any(x -> x <= 0.0 || !isfinite(x), eigC)
        return true
    end

    gap12 = abs(eigC[1] - eigC[2])
    gap23 = abs(eigC[2] - eigC[3])
    gap13 = abs(eigC[1] - eigC[3])

    scale = max(maximum(abs.(eigC)), 1.0)

    return min(gap12, gap23, gap13) < tol * scale
end

function Ψ(F_arg, Θ, p_field, mp::Ogden)

    c1, c2, m1, m2, K = mp.c1, mp.c2, mp.m1, mp.m2, mp.K

    if Θ <= 0.0 || !isfinite(Θ) || !isfinite(p_field)
        return 1.0e30
    end

    C = tdot(F_arg)
    Cmat = Matrix(C)

    eigC = eigvals(Symmetric(Cmat))

    if any(x -> x <= 0.0 || !isfinite(x), eigC)
        return 1.0e30
    end

    λ1 = sqrt(eigC[1])
    λ2 = sqrt(eigC[2])
    λ3 = sqrt(eigC[3])

    J = λ1 * λ2 * λ3

    if J <= 0.0 || !isfinite(J)
        return 1.0e30
    end

    Jm13 = J^(-1.0 / 3.0)

    λ̃1 = Jm13 * λ1
    λ̃2 = Jm13 * λ2
    λ̃3 = Jm13 * λ3

    if λ̃1 <= 0.0 || λ̃2 <= 0.0 || λ̃3 <= 0.0 ||
       !isfinite(λ̃1) || !isfinite(λ̃2) || !isfinite(λ̃3)

        return 1.0e30
    end

    Ψiso =
        c1 / m1^2 * (λ̃1^m1 + λ̃2^m1 + λ̃3^m1 - 3.0) +
        c2 / m2^2 * (λ̃1^m2 + λ̃2^m2 + λ̃3^m2 - 3.0)

    Ψvol = 0.5 * K * (log(Θ))^2

    Ψcoup = p_field * (J - Θ)

    Ψtotal = Ψiso + Ψvol + Ψcoup

    if !isfinite(Ψtotal)
        return 1.0e30
    end

    return Ψtotal
end

function numerical_P(F, Θ, p_field, mp::Ogden; εP=1.0e-6)

    Pmat = zeros(3, 3)

    for a in 1:3
        for b in 1:3
            Eab = basis_tensor_2(a, b)

            Ψplus = Ψ(F + εP * Eab, Θ, p_field, mp)
            Ψminus = Ψ(F - εP * Eab, Θ, p_field, mp)

            Pmat[a, b] = (Ψplus - Ψminus) / (2.0 * εP)
        end
    end

    return tensor2_from_matrix(Pmat)
end

function numerical_tangent_F(F, Θ, p_field, mp::Ogden; εA=1.0e-5, εP=1.0e-6)

    A = zeros(Float64, 3, 3, 3, 3)

    for a in 1:3
        for b in 1:3
            Eab = basis_tensor_2(a, b)

            Pplus = numerical_P(F + εA * Eab, Θ, p_field, mp; εP=εP)
            Pminus = numerical_P(F - εA * Eab, Θ, p_field, mp; εP=εP)

            dP = (Pplus - Pminus) / (2.0 * εA)

            for i in 1:3
                for j in 1:3
                    A[i, j, a, b] = dP[i, j]
                end
            end
        end
    end

    return tensor4_from_array(A)
end

function constitutive_driver(F_arg, Θ, p_field, mp::Ogden)

    use_numeric = has_repeated_eigs(F_arg)

    if !use_numeric
        try
            ∂²Ψ∂F²_ad, ∂Ψ∂F_ad =
                Tensors.hessian(y -> Ψ(y, Θ, p_field, mp), F_arg, :all)

            if finite_tensor(∂Ψ∂F_ad) && finite_tensor(∂²Ψ∂F²_ad)
                ∂Ψ∂F = ∂Ψ∂F_ad
                ∂²Ψ∂F² = ∂²Ψ∂F²_ad
            else
                use_numeric = true
            end
        catch
            use_numeric = true
        end
    end

    if use_numeric
        ∂Ψ∂F = numerical_P(F_arg, Θ, p_field, mp)
        ∂²Ψ∂F² = numerical_tangent_F(F_arg, Θ, p_field, mp)
    end

    K = mp.K

    if Θ <= 0.0 || !isfinite(Θ)
        ∂Ψ∂Θ = 1.0e30
        ∂²Ψ∂Θ² = 1.0e30
    else
        lnΘ = log(Θ)
        ∂Ψ∂Θ = K * lnΘ / Θ - p_field
        ∂²Ψ∂Θ² = K * (1.0 - lnΘ) / Θ^2
    end

    return ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂Θ, ∂²Ψ∂Θ²
end

# ================================================================================
# Simple gravity body force
# ================================================================================

function static_moving_frame_gravity_body_force(ρ, A, load_factor)
    return -load_factor * ρ * A
end

# ================================================================================
# Element assembly
# ================================================================================
function assemble_element!(Ke, fe, cell, cvu, cvp, cvΘ, mp, ue, pe, Θe, ρ, A, load_factor)

    reinit!(cvu, cell)
    reinit!(cvp, cell)
    reinit!(cvΘ, cell)

    fill!(Ke, 0.0)
    fill!(fe, 0.0)

    nu = getnbasefunctions(cvu)
    np = getnbasefunctions(cvp)
    nΘ = getnbasefunctions(cvΘ)

    b = static_moving_frame_gravity_body_force(ρ, A, load_factor)

    for qp in 1:getnquadpoints(cvu)

        dΩ = getdetJdV(cvu, qp)

        ∇u = function_gradient(cvu, qp, ue)
        F_qp = one(∇u) + ∇u

        p̂ = function_value(cvp, qp, pe)
        Θ̂ = function_value(cvΘ, qp, Θe)

        ∂Ψ∂F, ∂²Ψ∂F², ∂Ψ∂Θ, ∂²Ψ∂Θ² =
            constitutive_driver(F_qp, Θ̂, p̂, mp)

        if !finite_tensor(∂Ψ∂F) ||
           !finite_tensor(∂²Ψ∂F²) ||
           !isfinite(∂Ψ∂Θ) ||
           !isfinite(∂²Ψ∂Θ²)

            println("\nNaN/Inf from constitutive_driver")
            println("cell = ", cellid(cell))
            println("qp = ", qp)
            println("detF = ", det(F_qp))
            println("Θ = ", Θ̂)
            println("p = ", p̂)
            println("eigC = ", eigvals(Symmetric(Matrix(tdot(F_qp)))))
            error("Stopping: bad constitutive output")
        end

        Finv = inv(F_qp)
        J = det(F_qp)

        if !isfinite(J) || !finite_tensor(Finv)
            println("\nNaN/Inf from det(F) or inv(F)")
            println("cell = ", cellid(cell))
            println("qp = ", qp)
            println("F = ", F_qp)
            error("Stopping: bad F")
        end

        # ----------------------------------------------------------------------
        # Block 1: displacement equilibrium
        # ----------------------------------------------------------------------

        for i in 1:nu

            ∇Nᵢ = shape_gradient(cvu, qp, i)
            Nᵢ = shape_value(cvu, qp, i)

            # Internal residual
            fe[BlockIndex((1,), (i,))] += (∇Nᵢ ⊡ ∂Ψ∂F) * dΩ

            # External gravity body force
            fe[BlockIndex((1,), (i,))] -= (Nᵢ ⋅ b) * dΩ

            # Kuu
            for j in 1:nu
                ∇Nⱼ = shape_gradient(cvu, qp, j)

                Ke[BlockIndex((1, 1), (i, j))] +=
                    (∇Nᵢ ⊡ ∂²Ψ∂F² ⊡ ∇Nⱼ) * dΩ
            end

            # δJ_i
            δJ_i = J * tr(Finv ⋅ ∇Nᵢ)

            # Kup and Kpu
            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)

                Ke[BlockIndex((1, 2), (i, j))] += δJ_i * Nⱼᵖ * dΩ
                Ke[BlockIndex((2, 1), (j, i))] += δJ_i * Nⱼᵖ * dΩ
            end
        end

        # ----------------------------------------------------------------------
        # Block 2: pressure equation
        # ----------------------------------------------------------------------

        for i in 1:np

            Nᵢᵖ = shape_value(cvp, qp, i)

            fe[BlockIndex((2,), (i,))] += Nᵢᵖ * (J - Θ̂) * dΩ

            for j in 1:nΘ
                NⱼΘ = shape_value(cvΘ, qp, j)

                Ke[BlockIndex((2, 3), (i, j))] -= Nᵢᵖ * NⱼΘ * dΩ
            end
        end

        # ----------------------------------------------------------------------
        # Block 3: Θ equation
        # ----------------------------------------------------------------------

        for i in 1:nΘ

            NᵢΘ = shape_value(cvΘ, qp, i)

            fe[BlockIndex((3,), (i,))] += NᵢΘ * ∂Ψ∂Θ * dΩ

            for j in 1:nΘ
                NⱼΘ = shape_value(cvΘ, qp, j)

                Ke[BlockIndex((3, 3), (i, j))] +=
                    NᵢΘ * ∂²Ψ∂Θ² * NⱼΘ * dΩ
            end

            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)

                Ke[BlockIndex((3, 2), (i, j))] -= NᵢΘ * Nⱼᵖ * dΩ
            end
        end
    end

    return
end

# ================================================================================
# Global assembly
# ================================================================================

function assemble_global!( K::SparseMatrixCSC, f, cvu, cvp , cvΘ, dh, mp, w, ρ, A, load_factor)

    nu = getnbasefunctions(cvu)
    np = getnbasefunctions(cvp)
    nΘ = getnbasefunctions(cvΘ)
    ntot = nu + np + nΘ

    fe = BlockedArray(zeros(ntot), [nu, np, nΘ])
    ke = BlockedArray(zeros(ntot, ntot), [nu, np, nΘ], [nu, np, nΘ])

    assembler = start_assemble(K, f)

    for cell in CellIterator(dh)

        global_dofs = celldofs(cell)

        global_dofs_u = global_dofs[1:nu]
        global_dofs_p = global_dofs[(nu + 1):(nu + np)]
        global_dofs_Θ = global_dofs[(nu + np + 1):end]

        @assert length(global_dofs) == ntot

        ue = w[global_dofs_u]
        pe = w[global_dofs_p]
        Θe = w[global_dofs_Θ]

        assemble_element!( ke, fe, cell, cvu, cvp, cvΘ, mp, ue, pe , Θe, ρ, A, load_factor)
        assemble!(assembler, global_dofs, ke, fe)
    end

    return
end

# ================================================================================
# FE setup
# ================================================================================

function create_values(ipu, ipp, ipΘ)
    qr = QuadratureRule{RefTetrahedron}(2)

    cvu = CellValues(qr, ipu)
    cvp = CellValues(qr, ipp)
    cvΘ = CellValues(qr, ipΘ)

    return cvu, cvp, cvΘ
end

function create_dofhandler(grid, ipu, ipp, ipΘ)
    dh = DofHandler(grid)

    add!(dh, :u, ipu)
    add!(dh, :p, ipp)
    add!(dh, :Θ, ipΘ)

    close!(dh)

    return dh
end

# ================================================================================
# Solve
# ================================================================================

function solve(c1, c2, m1, m2, bulk, grid, load_time, numSteps, tissueDensity, gravityVector; NEWTON_TOL = 1.0e-8, NEWTON_MAXITER = 100)

    mp = Ogden(c1, c2, m1, m2, bulk)

    ipu = Lagrange{RefTetrahedron, 1}()^3
    ipp = DiscontinuousLagrange{RefTetrahedron, 0}()
    ipΘ = DiscontinuousLagrange{RefTetrahedron, 0}()

    dh = create_dofhandler(grid,  ipu, ipp, ipΘ)
    dbcs = create_bc(dh)
    cvu, cvp, cvΘ = create_values(ipu, ipp, ipΘ)

    nd = ndofs(dh)

    UT = Vector{Vector{Point{3, Float64}}}(undef, numSteps + 1)
    UT_mag = Vector{Vector{Float64}}(undef, numSteps + 1)
    ut_mag_max = zeros(Float64, numSteps + 1)

    un = zeros(nd)
    u = zeros(nd)
    Δu = zeros(nd)
    ΔΔu = zeros(nd)

    Tf = load_time
    Δt = Tf / numSteps

    # Initialize Θ DOFs to 1.0
    Θ_dof_range = dof_range(dh, :Θ)

    for cell in CellIterator(dh)
        gdofs = celldofs(cell)

        for d in gdofs[Θ_dof_range]
            un[d] = 1.0
        end
    end

    apply!(un, dbcs)
    u .= un

    # Store initial state
    u_nodes = vec(evaluate_at_grid_nodes(dh, u, :u))

    ux = getindex.(u_nodes, 1)
    uy = getindex.(u_nodes, 2)
    uz = getindex.(u_nodes, 3)

    disp_pts = [
        Point{3, Float64}((ux[j], uy[j], uz[j]))
        for j in eachindex(ux)
    ]

    UT[1] = disp_pts
    UT_mag[1] = norm.(disp_pts)
    ut_mag_max[1] = maximum(UT_mag[1])

    ρ = tissueDensity
    A = gravityVector

    for step in 1:numSteps

        t = step * Δt
        load_factor = step / numSteps

        println("\n=== Load step $step / $numSteps, t = $t, load_factor = $load_factor ===")

        Ferrite.update!(dbcs, t)

        Δu .= 0.0
        newton_itr = 0

        while true

            u .= un .+ Δu
            apply!(u, dbcs)

            Ksys = allocate_matrix(dh)
            g = zeros(nd)

            assemble_global!( Ksys, g, cvu, cvp, cvΘ, dh, mp, u, ρ, A, load_factor)

            normg = norm(g[Ferrite.free_dofs(dbcs)])

            println("  Newton iteration $newton_itr: residual = $normg")

            if normg < NEWTON_TOL
                println("  Converged in $newton_itr iterations")
                break
            elseif newton_itr ≥ NEWTON_MAXITER
                error("Newton failed to converge at t = $t")
            end

            free = Ferrite.free_dofs(dbcs)

            if any(!isfinite, Ksys.nzval) || any(!isfinite, g)
                println("Ksys finite = ", all(isfinite, Ksys.nzval))
                println("g finite = ", all(isfinite, g))
                error("Stopping: Ksys or g contains NaN/Inf")
            end

            ΔΔu .= 0.0

            ΔΔu[free] .= Ksys[free, free] \ g[free]

            α = 1.0
            Δu .-= α .* ΔΔu

            newton_itr += 1
        end

        un .= u

        u_nodes = vec(evaluate_at_grid_nodes(dh, u, :u))

        ux = getindex.(u_nodes, 1)
        uy = getindex.(u_nodes, 2)
        uz = getindex.(u_nodes, 3)

        disp_pts = [
            Point{3, Float64}((ux[j], uy[j], uz[j]))
            for j in eachindex(ux)
        ]

        UT[step + 1] = disp_pts
        UT_mag[step + 1] = norm.(disp_pts)
        ut_mag_max[step + 1] = maximum(UT_mag[step + 1])
    end

    return UT, UT_mag, ut_mag_max
end

UT, UT_mag, ut_mag_max = solve( c1, c2 , m1, m2, Kbulk, grid, load_time, numSteps, tissueDensity,  gravityVector)



# ================================================================================
# Visualization
# ================================================================================

numInc = length(UT)

scale = 1.0

VT = [
    V .+ scale .* UT[i]
    for i in 1:numInc
]

incRange = 0:1:(numInc - 1)

min_p = minp([minp(VT_i) for VT_i in VT])
max_p = maxp([maxp(VT_i) for VT_i in VT])

GLMakie.closeall()

fig_disp = Figure(size = (1000, 700))

stepStart = incRange[end]

ax3 = AxisGeom(
    fig_disp[1, 1],
    title = "Step: $stepStart",
    limits = (
        min_p[1],
        max_p[1],
        min_p[2],
        max_p[2],
        min_p[3],
        max_p[3],
    ),
)

hp = meshplot!(
    ax3,
    Fb,
    VT[stepStart + 1];
    strokewidth = 2,
    color = UT_mag[stepStart + 1],
    transparency = false,
    colormap = Reverse(:Spectral),
    colorrange = (0, maximum(ut_mag_max)),
)

Colorbar(
    fig_disp[1, 2],
    hp.plots[1],
    label = "Displacement magnitude [mm]",
)

hSlider = Slider(
    fig_disp[2, 1],
    range = incRange,
    startvalue = stepStart,
    linewidth = 30,
)

on(hSlider.value) do stepIndex
    i = stepIndex + 1

    hp[1] = GeometryBasics.Mesh(VT[i], Fb)
    hp.color = UT_mag[i]

    ax3.title = "Step: $stepIndex"
end

slidercontrol(hSlider, ax3)

display(GLMakie.Screen(), fig_disp)