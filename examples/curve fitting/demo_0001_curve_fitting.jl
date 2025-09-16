using Revise
using FerriteHyperelastic
a = [0	0
0.392330006	0.29642895
0.632225693	0.701048497
0.737471477	0.967418619
0.947946539	1.519877146
1.216990514	2.102009378
1.42757286	2.526299067
1.685063153	2.911229215
1.954354706	3.197587206
2.17095322	3.434562838
2.411022212	3.632140416]

inputDataType = "strain"
data_type     = "biaxial"
modelType     = "mooney-rivlin"
Sexp = a[:, 1]

if inputDataType == "lambda"
    strainExp = a[:, 2] .- 1
elseif inputDataType == "strain"
    strainExp = a[:, 2]
else
    error("Unknown inputDataType: $inputDataType")
end

mat_cons_solver = solver_constants_hyper(data_type,modelType,strainExp,Sexp)


