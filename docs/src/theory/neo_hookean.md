
## Compressible (unconstrained) neo-Hookean 

The strain energy for compressible neo-Hookean model is  
 

```math
\begin{equation}
\label{eq:1}
\Psi = \frac{\mu}{2}(I_1 - 3) - \mu \log J + \frac{\lambda}{2}(\log J)^2
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
Jacobian:

```math
J = \sqrt{\det(\mathbf{C})}
```

Material Parameters:

Young’s modulus: $E$

Poisson’s ratio: $\nu$

```math
\mu =  \frac{E}{2(1+\nu)}
```

```math
\lambda = \frac{\nu E}{(1+\nu)(1-2\nu)}
```

## Nearly-Incompressible (uncoupled deviatoric) neo-Hookean 

```math
\begin{equation}
\label{eq:2}
\Psi = C_1(\tilde{I}_1 - 3) + \frac{1}{2} K (\log J)^2
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

!!! note
    ```math
    \tilde{I}_1 = J^{-2/3}I_1
    ```
   