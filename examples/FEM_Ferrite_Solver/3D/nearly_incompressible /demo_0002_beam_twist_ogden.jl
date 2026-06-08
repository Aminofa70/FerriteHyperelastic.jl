using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics
using IterativeSolvers
using ComodoFerrite
using ComodoFerrite.Ferrite
using SparseArrays
using BlockArrays
using Tensors
using LinearAlgebra
using DelimitedFiles

## Geometry and Mesh Generation
## Hex_20 mesh type 
function create_grid(Lx, Ly, Lz, nx, ny, nz)
    left = Ferrite.Vec(0.0, 0.0, 0.0)
    right = Ferrite.Vec(Lx, Ly, Lz)
    grid = generate_grid(Hexahedron, (nx, ny, nz), left, right)
    return grid
end

### call the grid function and generate the mesh
Lx, Ly, Lz = 10.0, 40.0, 10.0
nx, ny, nz = 5, 20, 5
grid = create_grid(Lx, Ly, Lz, nx, ny, nz)

## define the material parameters
c1 = 1e-3 
m1 = 8
c2 = c1
m2 = -m1
k_factor = 1e2  #Bulk modulus factor 
k = c1*k_factor #Bulk modulus

##  define the load factor and number of step
load_factor = 1.0      # 1.0 = full 360 degree twist
numSteps = 25

# Usual 2D rotation about the origin in the x-z plane:
#
# [x_new]   [ cos(θ)   sin(θ)] [x]
# [z_new] = [-sin(θ)   cos(θ)] [z]
#
# Therefore:
# x_new = x*cos(θ) + z*sin(θ)
# z_new = -x*sin(θ) + z*cos(θ)


# the rotation is about the center (cx,cz​)
# clockwise rotation by an angle θ
# Rotation about the center (cx, cz) in the x-z plane:
#
# [x_new]   [cx] + [ cos(θ)   sin(θ)] [x - cx]
# [z_new] = [cz] + [-sin(θ)   cos(θ)] [z - cz]
#
# Therefore:
# x_new = cx + (x - cx)*cos(θ) + (z - cz)*sin(θ)
# z_new = cz - (x - cx)*sin(θ) + (z - cz)*cos(θ)

function rotation_front(X, t; Lx = 10.0, Lz = 10.0)
    θ = 2π * t  # full 360° twist, like the GIBBON demo
    cx, cz = Lx / 2, Lz / 2
    x, y, z = X
    x_new = cx + (x - cx) * cos(θ) + (z - cz) * sin(θ)
    z_new = cz - (x - cx) * sin(θ) + (z - cz) * cos(θ)
    return Ferrite.Vec{3}((x_new - x, 0.0, z_new - z))
end

## create Dirichlet boundary conditions
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "back"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3]))
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "front"), (x, t) -> rotation_front(x, t), [1, 2, 3]))
    Ferrite.close!(ch)
    Ferrite.update!(ch, 0.0)
    return ch
end

## material modeling (Ogden)
## Derive the gradient (pressure) and Hessian (Tangent stiffness)
struct Ogden
    c1::Float64
    c2::Float64
    m1::Float64
    m2::Float64
    K::Float64
end

# ══════════════════════════════════════════════════════════════
# HYBRID OGDEN CONSTITUTIVE DRIVER
# Fast path: Tensors.hessian when principal stretches are distinct.
# Safe path: numerical stress/tangent when eigenvalues of C are repeated
#            or when AD Hessian returns NaN/Inf.
#
# Reason:
# AD Hessian through eigvals is not reliable at repeated eigenvalues,
# e.g. C = I -> eigC = [1, 1, 1].
# ══════════════════════════════════════════════════════════════

# Converts a 3×3 real-valued matrix into a second-order tensor in 3D.
# The input matrix A is copied into a Tensor{2, 3} using Float64 values.
function tensor2_from_matrix(A::AbstractMatrix{<:Real})
    return Tensor{2, 3}(Float64[
        A[1, 1] A[1, 2] A[1, 3]
        A[2, 1] A[2, 2] A[2, 3]
        A[3, 1] A[3, 2] A[3, 3]
    ])
end

# Converts a 3×3×3×3 Float64 array into a fourth-order tensor in 3D.
# The input array A is used directly to construct a Tensor{4, 3}.
function tensor4_from_array(A::Array{Float64, 4})
    return Tensor{4, 3}(A)
end

# Creates a second-order basis tensor e_a ⊗ e_b in 3D.
# The returned tensor has value 1.0 at index (a, b) and 0.0 everywhere else.
function basis_tensor_2(a::Int, b::Int)
    A = zeros(3, 3)
    A[a, b] = 1.0
    return tensor2_from_matrix(A)
end

# Checks whether all components of a tensor-like object are finite.
# Returns true only if every value in x is finite, meaning no NaN or Inf values are present.
function finite_tensor(x)
    return all(isfinite, collect(x))
end

# Detect when the spectral AD Hessian is dangerous.
# If two eigenvalues are equal/very close, the ordered eigenvalue map
# is not smooth enough for a reliable Hessian.
function has_repeated_eigs(F; tol = 1.0e-8)
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

    # Check before sqrt.
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

# Numerical first Piola stress:
# P_ab = ∂Ψ / ∂F_ab
function numerical_P(F, Θ, p_field, mp::Ogden; εP = 1.0e-6)

    Pmat = zeros(3, 3)

    for a in 1:3
        for b in 1:3
            Eab = basis_tensor_2(a, b)

            Ψplus  = Ψ(F + εP * Eab, Θ, p_field, mp)
            Ψminus = Ψ(F - εP * Eab, Θ, p_field, mp)

            Pmat[a, b] = (Ψplus - Ψminus) / (2.0 * εP)
        end
    end

    return tensor2_from_matrix(Pmat)
end

# Numerical material tangent:
# A_ijab = ∂P_ij / ∂F_ab
function numerical_tangent_F(F, Θ, p_field, mp::Ogden; εA = 1.0e-5, εP = 1.0e-6)

    A = zeros(Float64, 3, 3, 3, 3)

    for a in 1:3
        for b in 1:3
            Eab = basis_tensor_2(a, b)

            Pplus  = numerical_P(F + εA * Eab, Θ, p_field, mp; εP = εP)
            Pminus = numerical_P(F - εA * Eab, Θ, p_field, mp; εP = εP)

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

function assemble_element!(Ke, fe, cell, cvu, cvp, cvΘ,mp, ue, pe, Θe)

    reinit!(cvu, cell)
    reinit!(cvp, cell)
    reinit!(cvΘ, cell)

    fill!(Ke, 0.0)
    fill!(fe, 0.0)

    nu = getnbasefunctions(cvu)
    np = getnbasefunctions(cvp)
    nΘ = getnbasefunctions(cvΘ)

    for qp in 1:getnquadpoints(cvu)
        dΩ = getdetJdV(cvu, qp)

        # --- Current state at quadrature point ---
        ∇u = function_gradient(cvu, qp, ue)
        F_qp = one(∇u) + ∇u
        p̂ = function_value(cvp, qp, pe)
        Θ̂ = function_value(cvΘ, qp, Θe)

        # --- Constitutive evaluation ---
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
        # ════════════════════════════════════════════════════
        # Block 1 — u-equation (equilibrium)
        #   Rᵤ[i] = ∫ ∇Nᵢ : P dV,  P = ∂Ψ/∂F
        # ════════════════════════════════════════════════════
        for i in 1:nu
            ∇Nᵢ = shape_gradient(cvu, qp, i)

            # Residual
            fe[BlockIndex((1,), (i,))] += (∇Nᵢ ⊡ ∂Ψ∂F) * dΩ

            # Kuu — material + geometric stiffness
            for j in 1:nu
                ∇Nⱼ = shape_gradient(cvu, qp, j)
                Ke[BlockIndex((1, 1), (i, j))] += (∇Nᵢ ⊡ ∂²Ψ∂F² ⊡ ∇Nⱼ) * dΩ
            end

            # δJ_i = J F⁻ᵀ : ∇Nᵢ  (variation of det F)
            δJ_i = J * tr(Finv ⋅ ∇Nᵢ)

            # Kup & Kpu — symmetric coupling
            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)
                Ke[BlockIndex((1, 2), (i, j))] += δJ_i * Nⱼᵖ * dΩ
                Ke[BlockIndex((2, 1), (j, i))] += δJ_i * Nⱼᵖ * dΩ
            end
        end

        # ════════════════════════════════════════════════════
        # Block 2 — p-equation (constraint J = Θ)
        #   Rₚ[i] = ∫ Nᵢᵖ (J − Θ) dV
        # ════════════════════════════════════════════════════
        for i in 1:np
            Nᵢᵖ = shape_value(cvp, qp, i)

            # Residual
            fe[BlockIndex((2,), (i,))] += Nᵢᵖ * (J - Θ̂) * dΩ

            # KpΘ — off-diagonal coupling
            for j in 1:nΘ
                NⱼΘ = shape_value(cvΘ, qp, j)
                Ke[BlockIndex((2, 3), (i, j))] -= Nᵢᵖ * NⱼΘ * dΩ
            end
        end

        # ════════════════════════════════════════════════════
        # Block 3 — Θ-equation (volumetric constitutive)
        #   RΘ[i] = ∫ NᵢΘ (U'(Θ) − p) dV
        #   KΘΘ uses U''(Θ) = K(1 − ln Θ)/Θ²
        # ════════════════════════════════════════════════════
        #
        for i in 1:nΘ
            NᵢΘ = shape_value(cvΘ, qp, i)

            # Residual
            fe[BlockIndex((3,), (i,))] += NᵢΘ * ∂Ψ∂Θ * dΩ

            # KΘΘ — volumetric tangent
            for j in 1:nΘ
                NⱼΘ = shape_value(cvΘ, qp, j)
                Ke[BlockIndex((3, 3), (i, j))] += NᵢΘ * ∂²Ψ∂Θ² * NⱼΘ * dΩ
            end

            # KΘp — (= KpΘᵀ)
            for j in 1:np
                Nⱼᵖ = shape_value(cvp, qp, j)
                Ke[BlockIndex((3, 2), (i, j))] -= NᵢΘ * Nⱼᵖ * dΩ
            end
        end
    end
    return
end

# ══════════════════════════════════════════════════════════════
# 7. GLOBAL ASSEMBLY
# ══════════════════════════════════════════════════════════════

function assemble_global!(K::SparseMatrixCSC, f, cvu, cvp, cvΘ, dh, mp, w)

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

        assemble_element!(ke, fe, cell, cvu, cvp, cvΘ, mp, ue, pe, Θe)

        assemble!(assembler, global_dofs, ke, fe)
    end
    return
end

# ══════════════════════════════════════════════════════════════
# 8. FE SETUP HELPERS
# ══════════════════════════════════════════════════════════════

function create_values(ipu, ipp, ipΘ)
    qr = QuadratureRule{RefHexahedron}(2)
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


function get_front_nodes(grid, Ly)
    nodes = getnodes(grid)
    return [i for i in eachindex(nodes) if nodes[i].x[2] ≈ Ly]
end

function compute_front_torque_y(dh, grid, residual, front_nodes, Lx, Lz)
    f_nodes = vec(evaluate_at_grid_nodes(dh, residual, :u))
    nodes = getnodes(grid)

    cx = Lx / 2
    cz = Lz / 2

    My = 0.0

    for node in front_nodes
        X = nodes[node].x

        rx = X[1] - cx
        rz = X[3] - cz

        Fx = f_nodes[node][1]
        Fz = f_nodes[node][3]

        My += rx * Fz - rz * Fx
    end

    return My
end

#Standard Newton iteration for the three-field principle.
function solve( c1, c2, m1, m2, bulk, grid, load_factor, numSteps; NEWTON_TOL = 1.0e-8, NEWTON_MAXITER = 100)

    # --- Material ---
    mp = Ogden(c1, c2, m1, m2, bulk)


    # --- Interpolations: Q1/P0/P0 ---
    ipu = Lagrange{RefHexahedron, 1}()^3
    ipp = DiscontinuousLagrange{RefHexahedron, 0}()
    ipΘ = DiscontinuousLagrange{RefHexahedron, 0}()

    # --- FE setup ---
    dh = create_dofhandler(grid, ipu, ipp, ipΘ)
    dbcs = create_bc(dh)
    cvu, cvp, cvΘ = create_values(ipu, ipp, ipΘ)

    dbcs = create_bc(dh)
    nd = ndofs(dh)

    # --- Front nodes for torque ---
    front_nodes = get_front_nodes(grid, Ly)

    # --- Storage for visualization ---
    UT = Vector{Vector{Point{3, Float64}}}(undef, numSteps + 1)
    UT_mag = Vector{Vector{Float64}}(undef, numSteps + 1)
    ut_mag_max = zeros(Float64, numSteps + 1)

    # --- Storage for torque plot ---
    My_time = zeros(Float64, numSteps + 1)
    My_front = zeros(Float64, numSteps + 1)

    # --- Solution vectors ---
    un = zeros(nd)
    u = zeros(nd)
    Δu = zeros(nd)
    ΔΔu = zeros(nd)

    # --- Time/load parameters ---
    Tf = load_factor
    Δt = Tf / numSteps

    # --- Initialize Θ DOFs to 1.0 ---
    Θ_dof_range = dof_range(dh, :Θ)

    for cell in CellIterator(dh)
        gdofs = celldofs(cell)
        for d in gdofs[Θ_dof_range]
            un[d] = 1.0
        end
    end

    apply!(un, dbcs)
    u .= un

    # --- Store initial state ---
    u_nodes = vec(evaluate_at_grid_nodes(dh, u, :u))
    ux = getindex.(u_nodes, 1)
    uy = getindex.(u_nodes, 2)
    uz = getindex.(u_nodes, 3)

    disp_pts = [
        Point{3, Float64}([ux[j], uy[j], uz[j]])
        for j in eachindex(ux)
    ]

    UT[1] = disp_pts
    UT_mag[1] = norm.(disp_pts)
    ut_mag_max[1] = maximum(UT_mag[1])

    My_time[1] = 0.0
    My_front[1] = 0.0

    # ══════════════════════════════════════════════════════════
    # TIME STEPPING + NEWTON
    # ══════════════════════════════════════════════════════════
    for step in 1:numSteps

        t = step * Δt
        My_time[step + 1] = t

        println("\n=== Time step $step / $numSteps,  t = $t ===")

        Ferrite.update!(dbcs, t)

        Δu .= 0.0
        newton_itr = 0

        while true

            u .= un .+ Δu
            apply!(u, dbcs)

            # Fresh system every Newton iteration
            Ksys = allocate_matrix(dh)
            g = zeros(nd)

            assemble_global!(Ksys, g,cvu, cvp, cvΘ,dh, mp, u)

            normg = norm(g[Ferrite.free_dofs(dbcs)])

            if normg < NEWTON_TOL
                println("  Converged in $newton_itr iterations  (‖g‖ = $normg)")
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

            # Optional damping. Keep 1.0 for full Newton.
            α = 1.0
            Δu .-= α .* ΔΔu

            newton_itr += 1
        end

        # Accept converged solution
        un .= u

        # ══════════════════════════════════════════════════════
        # Reaction torque on front face
        # Use fresh residual, before applying zero constraints
        # ══════════════════════════════════════════════════════
        K_react = allocate_matrix(dh)
        g_react = zeros(nd)

        assemble_global!(K_react,g_react,cvu, cvp, cvΘ,dh, mp, u)

        My_front[step + 1] = compute_front_torque_y( dh, grid, g_react,front_nodes, Lx, Lz)

        println("  My_front = ", My_front[step + 1])

        # --- Post-processing displacement ---
        u_nodes = vec(evaluate_at_grid_nodes(dh, u, :u))
        ux = getindex.(u_nodes, 1)
        uy = getindex.(u_nodes, 2)
        uz = getindex.(u_nodes, 3)

        disp_pts = [
            Point{3, Float64}([ux[j], uy[j], uz[j]])
            for j in eachindex(ux)
        ]

        UT[step + 1] = disp_pts
        UT_mag[step + 1] = norm.(disp_pts)
        ut_mag_max[step + 1] = maximum(UT_mag[step + 1])
    end

    return UT, UT_mag, ut_mag_max, My_time, My_front
end

load_factor = 1.0
numSteps = 25

UT, UT_mag, ut_mag_max, time_curve, My = solve( c1, c2, m1, m2, k,grid, load_factor, numSteps)

E , V, F, Fb, CFb_type   = FerriteToComodo(grid)



numInc = length(UT)

# Create displaced mesh per step
scale = 1.0
VT = [V .+ scale .* UT[i] for i in 1:(numSteps + 1)]
incRange = 0:1:(numInc - 1)

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

#######
# Visualization
GLMakie.closeall()

fig = Figure(size=(1200,1000))
stepStart = incRange[end]
ax1 = AxisGeom(fig[1, 1], title = "Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp1 = meshplot!(ax1, Fb, VT[end]; strokewidth=2, color=UT_mag[end], transparency=false, colormap = Reverse(:Spectral),colorrange=(0,maximum(ut_mag_max)))
Colorbar(fig[1, 2],hp1.plots[1],label = "Displacement magnitude [mm]") 

ax2 = Axis(fig[1, 3], title = "Step: $stepStart", xlabel="Time [s]", ylabel="Reaction torque [Nmm]")
lines!(ax2, time_curve, My, color=:red, linewidth=3)
hp2 = scatter!(ax2, Point{2,Float64}(time_curve[stepStart+1], My[stepStart+1]), markersize=15, color=:red)

hSlider = Slider(fig[2, :], range = incRange, startvalue = stepStart,linewidth=30)
on(hSlider.value) do stepIndex 
    hp1[1] = GeometryBasics.Mesh(VT[stepIndex+1],Fb)
    hp1.color = UT_mag[stepIndex+1]
    
    hp2[1] = Point{2,Float64}(time_curve[stepIndex+1], My[stepIndex+1])

    ax1.title = "Step: $stepIndex"
    ax2.title = "Step: $stepIndex"
end

slidercontrol(hSlider,ax1)

screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "FerriteHyperelastic example")