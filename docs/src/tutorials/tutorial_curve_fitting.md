# This tutorial is for curve fitting of hyperelastic models 
## We find the marterial parameters in hyperelastic models


For now the package supports neo-Hookean, Mooney-Rivlin, Ogden, and Yeoh models. 

The package also supports stability. Some experimental data are available in [HyperData.jl](https://github.com/Aminofa70/HyperData.jl) repository.  

The formulation for these models are found at [Abaqus tutorials](https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/stm/default.htm?startat=ch04s06ath123.html) . 

The experimental data can be stain or stretch for ```x ``` axis and the nominal stress for ```y ```axis.

Assume that we have experimental data in a matrix named ```data ```. 
Then we have the following code

```
inputDataType = "strain"
data_type     = "biaxial"
modelType     = "ogden"
Sexp = data[:, 1]

if inputDataType == "lambda"
    strainExp = data[:, 2] .- 1
elseif inputDataType == "strain"
    strainExp = data[:, 2]
else
    error("Unknown inputDataType: $inputDataType")
end

mat_cons_solver = solver_constants_hyper(data_type, modelType,strainExp, Sexp)

```

``` inputDataType ``` can be ```"lambda"``` or ```"strain"```.

Also ```  data_type ``` can be biaixla, uniaxial, and pure shear.

```modelType``` can be ogden, neo-hookean, yeoh, and mooney-rivlin. It is noted that we use lowercase letter for all.

