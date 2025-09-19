using Revise       # automatically reload changed code
using FerriteHyperelastic
using GLMakie
using HyperData

# Read experimental data (uniaxial test)
λ_exp, P_exp = read_data!(Treloar_1944, uniaxial)

# Input parameters
inputDataType = "lambda"
data_type = "uniaxial"
modelType = "yeoh"    # Changed to yeoh

Sexp = P_exp
strainExp = λ_exp .- 1

# Get material constants by fitting (assumed function)
mat_cons_solver = solver_constants_hyper(data_type, modelType, strainExp, Sexp)
C1, C2, C3 = mat_cons_solver[1], mat_cons_solver[2], mat_cons_solver[3]

# Generate strain values for smooth curve
ϵ = collect(LinRange(strainExp[1], strainExp[end], 100))

# Calculate stretches for incompressible uniaxial tension
λ1 = @. ϵ + 1
λ = λ1
λ2 = @. 1 / sqrt(λ)
λ3 = λ2   

# Calculate first invariant I1 for each λ
I1 = @. λ^2 + 2 / λ

# Yeoh nominal stress formula for uniaxial tension:
P_model = @. 2 * (λ - 1 / λ^2) * (C1 + 2 * C2 * (I1 - 3) + 3 * C3 * (I1 - 3)^2)

# Plot results
GLMakie.closeall()

fig = Figure(size=(800, 600), fontsize=26)
ax = Axis(fig[1, 1], xlabel= L"\mathscr{ε}", ylabel=L"P",  xgridvisible=false, ygridvisible=false)

lines!(ax, ϵ, P_model, color = :black, label="Fit, Yeoh")
scatter!(ax, strainExp, P_exp, marker=:circle, color=:red, label="Uniaxial experiment")

axislegend(ax, position=:lt, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(fig)
