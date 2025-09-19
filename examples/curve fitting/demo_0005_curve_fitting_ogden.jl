using Revise          # auto reload
using FerriteHyperelastic
using GLMakie
using HyperData

# Read experimental data (uniaxial test)
λ_exp, P_exp = read_data!(Treloar_1944, uniaxial)

# Input parameters
inputDataType = "lambda"
data_type = "uniaxial"
modelType = "ogden"   # <-- switched to Ogden

Sexp = P_exp
strainExp = λ_exp .- 1

# Fit material constants
# Expecting: [μ1, α1, μ2, α2, μ3, α3]
mat_cons_solver = solver_constants_hyper(data_type, modelType, strainExp, Sexp)
μ1, α1, μ2, α2, μ3, α3 = mat_cons_solver

# Smooth strain curve
ϵ = collect(LinRange(strainExp[1], strainExp[end], 200))

# Corresponding stretches
λ = @. ϵ + 1

# Ogden nominal stress (uniaxial, incompressible)
P_model = @. (2*μ1/α1)*(λ^(α1-1) - λ^(-1 - α1/2)) +
            (2*μ2/α2)*(λ^(α2-1) - λ^(-1 - α2/2)) +
            (2*μ3/α3)*(λ^(α3-1) - λ^(-1 - α3/2))

# Plot results
GLMakie.closeall()
fig = Figure(size=(800, 600), fontsize=26)
ax = Axis(fig[1, 1], xlabel=L"\varepsilon", ylabel=L"P", 
          xgridvisible=false, ygridvisible=false)

lines!(ax, ϵ, P_model, color=:black, label="Fit, Ogden (3-term)")
scatter!(ax, strainExp, P_exp, marker=:circle, color=:red, label="Uniaxial experiment")

axislegend(ax, position=:lt, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(fig)
