using Revise          # to automatically reload changed code
using FerriteHyperelastic
using GLMakie
using HyperData
# Read experimental data (uniaxial test)
λ_exp, P_exp = read_data!(Treloar_1944, equibiaxial)
# Input parameters
inputDataType = "lambda"
data_type = "biaxial"
modelType = "neo-hookean"
Sexp = P_exp
strainExp = λ_exp .- 1
# Get material constants 
mat_cons_solver = solver_constants_hyper(data_type, modelType, strainExp, Sexp)
C1 = mat_cons_solver[1]
# Generate strain values for smooth curve
ϵ = collect(LinRange(strainExp[1], strainExp[end], 100))
# Calculate stretches for incompressible uniaxial tension
λ1 = @. ϵ + 1
λ = λ1
λ2 = @. 1 / sqrt(λ)
λ3 = λ2   
# neo-Hookean stress (uniaxial nominal stress)
P_model = @. 2 * C1 * (λ - λ^-5)
# Plot results
GLMakie.closeall()
fig = Figure(size=(800, 600), fontsize=26)
ax = Axis(fig[1, 1], xlabel= L"\mathscr{ε}", ylabel=L"P",  xgridvisible=false, ygridvisible=false)
scatter!(ax, strainExp, P_exp, marker=:diamond, color=:red, label="Equibiaxial experiment", markersize = 16)
lines!(ax, ϵ, P_model, color = :black, label="Fit, neo-Hookean", linewidth = 3)

#############################################################################
# Read experimental data (uniaxial test)
λ_exp, P_exp = read_data!(Treloar_1944, equibiaxial)

# Input parameters
inputDataType = "lambda"
data_type = "biaxial"
modelType = "mooney-rivlin"    # Changed to mooney-rivlin

Sexp = P_exp
strainExp = λ_exp .- 1

# Get material constants by fitting (assumed function)
mat_cons_solver = solver_constants_hyper(data_type, modelType, strainExp, Sexp)
C1, C2 = mat_cons_solver[1], mat_cons_solver[2]

# Generate strain values for smooth curve
ϵ = collect(LinRange(strainExp[1], strainExp[end], 100))

# Calculate stretches for incompressible uniaxial tension
λ1 = @. ϵ + 1
λ = λ1
λ2 = @. 1 / sqrt(λ)
λ3 = λ2   
# Mooney-Rivlin nominal stress formula for uniaxial tension:
P_model = @. 2 * (C1 * (λ - 1/λ^5) + C2 * (λ^3 - 1/λ^3))
lines!(ax, ϵ, P_model, color = :blue, label="Fit, Mooney-Rivlin", linestyle = :dashdotdot, linewidth = 3)
#############################################################################
# Read experimental data (uniaxial test)
λ_exp, P_exp = read_data!(Treloar_1944, equibiaxial)

# Input parameters
inputDataType = "lambda"
data_type = "biaxial"
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

I1 =@. 2*λ^2 + λ^-4          # First invariant
P_model =@. 2 * (λ - λ^-5) * (C1 + 2*C2*(I1 - 3) + 3*C3*(I1 - 3)^2)

lines!(ax, ϵ, P_model, color = :green, label="Fit, Yeoh", linestyle = :dot, linewidth = 3)
#############################################################################
# Read experimental data (uniaxial test)
λ_exp, P_exp = read_data!(Treloar_1944, equibiaxial)

# Input parameters
inputDataType = "lambda"
data_type = "biaxial"
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
P_model = @. (2*μ1/α1)*(λ^(α1-1) - λ^(-2*α1-1)) +
            (2*μ2/α2)*(λ^(α2-1) - λ^(-2*α2-1)) +
            (2*μ3/α3)*(λ^(α3-1) - λ^(-2*α3-1))



lines!(ax, ϵ, P_model, color=:purple, label="Fit, Ogden (3-term)")
axislegend(ax, position=:lt, backgroundcolor=(:white, 0.7), framecolor=:gray,  linewidth = 3)

display(fig)
save("equibiaxial.png", fig) ## for save the plot