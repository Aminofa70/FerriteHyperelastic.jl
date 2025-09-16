using Revise
using FerriteHyperelastic
using GLMakie
using Ferrite
using Tensors
using HyperData
using LinearAlgebra
λ_exp, P_exp = read_data!(Treloar_1944, uniaxial)
inputDataType = "lambda"
data_type = "uniaxial"
modelType = "mooney-rivlin"

Sexp = P_exp
strainExp = λ_exp .- 1

# Get material constants by fitting (assumed function)
mat_cons_solver = solver_constants_hyper(data_type, modelType, strainExp, Sexp)

C1 = mat_cons_solver[1]
C2 = mat_cons_solver[2]

# Generate strain values for smooth curve
ϵ = collect(LinRange(strainExp[1], strainExp[end], 100))

# Calculate stretches for incompressible uniaxial tension
λ1 = @. ϵ + 1
λ = λ1
λ2 = @. 1 / sqrt(λ)
λ3 = λ2   # usually λ3 = λ2 for incompressibility (plane stress/strain assumption)

# Neo-Hookean stress (uniaxial nominal stress)
P_model = @. (2*C1 + ((2*C2)/λ))* (λ - 1 / λ^2)

# Plot results
GLMakie.closeall()

fig = Figure(size=(800, 600), fontsize=26)
ax = Axis(fig[1, 1], xlabel= L"\mathscr{ε}", ylabel=L"P",  xgridvisible=false, ygridvisible=false)

lines!(ax, ϵ, P_model, color = :black, label="Fit, Mooney-Rivlin")
scatter!(ax, strainExp, P_exp, marker=:circle, color=:red, label="Uniaxial experiment")

axislegend(ax, position=:lt, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(fig)


