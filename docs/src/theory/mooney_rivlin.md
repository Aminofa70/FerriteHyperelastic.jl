
## Comressible (unconstrained/coupled) Mooney-Rivlin

The strain energy for compressible Mooney-Rivlin model is  

```math
\begin{equation}
\label{eq:1}
\Psi = c_1(I_1 - 3) +c_2(I_2 -3) - 2(c_1+2c_2) \log J + \frac{K}{2}(\log J)^2
\end{equation}
```

Right Cauchy–Green tensor:
```math
\mathbf{C} = \mathbf{F}^\top \mathbf{F}
```
First invariant:
```math
I_1 = \mathrm{tr}(\mathbf{C})
```
Second invariant:
```math
I_2 = \tfrac{1}{2} \left[ (\mathrm{tr}\,\mathbf{C})^2 - \mathrm{tr}(\mathbf{C}^2) \right]
```

Jacobian:
```math
J = \sqrt{\det(\mathbf{C})}
```

## Nearly-Incompressible (uncoupled deviatoric) Mooney-Rivlin

```math
\begin{equation}
\label{eq:2}
\Psi = C_1(\tilde{I}_1 - 3) + C_2(\tilde{I}_2 - 3) +  \frac{1}{2} K (\log J)^2
\end{equation}
```

The determinant:

```math
J = \det( \mathbf{F})
```

The deformation gradient:

```math
\mathbf{\tilde{F}} = J^{-1/3} \mathbf{F}
``` 

The deviatoric part of the right Cauchy-Green deformation tensor:
```math
\mathbf{\tilde{C}} = \mathbf{\tilde{F}}^\top \mathbf{\tilde{F}}
```
The first invariant of the deviatoric part of the right Cauchy-Green deformation tensor:

```math
\tilde{I}_1 = \mathrm{tr}(\mathbf{\tilde{C}})
```

The second invariant of the deviatoric part of the right Cauchy-Green deformation tensor:

```math
\tilde{I}_2 = \tfrac{1}{2}\left[(\mathrm{tr}\,\tilde{\mathbf{C}})^2 - \mathrm{tr}(\tilde{\mathbf{C}}^2)\right]
```

!!! note
    ```math
    \tilde{I}_1 = J^{-2/3}I_1
    ```
    ```math
    \tilde{I}_2 = J^{-4/3}I_2
    ```