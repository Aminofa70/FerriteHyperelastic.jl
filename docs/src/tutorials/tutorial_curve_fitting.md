# This tutorial is for curve fitting in hyperelastic models.

!!! hint 
    The package supports 
    * neo-Hookean model
    * Mooney-Rivlin model
    * Ogden model
    * Yeoh models model

!!! note
    The package also supports *the Drucker stability criterion*.
    
In order to perform curve fitting, experimental data are required. 
Some experimental data are available in [HyperData.jl](https://github.com/Aminofa70/HyperData.jl) repository. 

The experimental data can be *Strain*  or  *Stretch* for ```x ``` axis and the nominal stress for ```y ```axis.


!!! note "Steps TODO Curve Fitting"
    * Load Package
    * Read Experimental Data
    * Curve Fitting Inputs
    * Call curve fitting function
    * Check the accuracy of the obtained material constant(s).

### Load Package

Load the required packages for curve fitting as
```
using Revise          # to automatically reload changed code
using FerriteHyperelastic
using GLMakie
using HyperData
```

### Read Experimental Data

Read the experimental data as
```
λ_exp, P_exp = read_data!(Treloar_1944, uniaxial)
```

### Curve Fitting Inputs

Set inputs of curve fitting function 

!!! warning
    The input MUST be Strain and the Nominal Stress, for this reasion if the data is different 
    use must convert them to starin and nominal stress

```
# Input parameters
inputDataType = "lambda"
data_type = "uniaxial"
modelType = "neo-hookean"

Sexp = P_exp
strainExp = λ_exp .- 1

```
!!! note 
    Please use lowercase letter for all model types (this example is "neo-hookean")

### Call curve fitting function
The curve-fitting function is now called. For neo-Hookean, the material constant is C1.
```
# Get material constants by fitting (assumed function)
mat_cons_solver = solver_constants_hyper(data_type, modelType, strainExp, Sexp)
C1 = mat_cons_solver[1]

```

### Check the accuracy of the obtained material constant(s)

```
# Generate strain values for smooth curve
ϵ = collect(LinRange(strainExp[1], strainExp[end], 100))

# Calculate stretches for incompressible uniaxial tension
λ1 = @. ϵ + 1
λ = λ1
λ2 = @. 1 / sqrt(λ)
λ3 = λ2   

# neo-Hookean stress (uniaxial nominal stress)
P_model = @. 2 * C1 * (λ - 1 / λ^2)

# Plot results
GLMakie.closeall()

fig = Figure(size=(800, 600), fontsize=26)
ax = Axis(fig[1, 1], xlabel= L"\mathscr{ε}", ylabel=L"P",  xgridvisible=false, ygridvisible=false)

lines!(ax, ϵ, P_model, color = :black, label="Fit, neo-Hookean")
scatter!(ax, strainExp, P_exp, marker=:circle, color=:red, label="Uniaxial experiment")

axislegend(ax, position=:lt, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(fig)


```

## Plain program
The version of the code without comments.

```
using Revise
using FerriteHyperelastic
using GLMakie
using HyperData


λ_exp, P_exp = read_data!(Treloar_1944, uniaxial)


inputDataType = "lambda"
data_type = "uniaxial"
modelType = "neo-hookean"

Sexp = P_exp
strainExp = λ_exp .- 1


mat_cons_solver = solver_constants_hyper(data_type, modelType, strainExp, Sexp)
C1 = mat_cons_solver[1]


ϵ = collect(LinRange(strainExp[1], strainExp[end], 100))


λ1 = @. ϵ + 1
λ = λ1
λ2 = @. 1 / sqrt(λ)
λ3 = λ2   


P_model = @. 2 * C1 * (λ - 1 / λ^2)


GLMakie.closeall()

fig = Figure(size=(800, 600), fontsize=26)
ax = Axis(fig[1, 1], xlabel= L"\mathscr{ε}", ylabel=L"P",  xgridvisible=false, ygridvisible=false)

lines!(ax, ϵ, P_model, color = :black, label="Fit, neo-Hookean")
scatter!(ax, strainExp, P_exp, marker=:circle, color=:red, label="Uniaxial experiment")

axislegend(ax, position=:lt, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(fig)
```
