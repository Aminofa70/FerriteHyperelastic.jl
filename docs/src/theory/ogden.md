
## Compressible (unconstrained) Ogden

The strain energy for compressible Ogden model is  

```math
\begin{equation}
\label{eq:1}
\Psi = \frac{1}{2}c_p (J-1)^2+ \sum_{i=1}^{N} \frac{c_i}{m_i} (\lambda_1^{m_i}+ \lambda_2^{m_i}+ \lambda_3^{m_i} -3  -m_i \log J)
\end{equation}
```

Right Cauchy–Green tensor:
```math
\mathbf{C} = \mathbf{F}^\top \mathbf{F}
```
The eigenvalues of right Cauchy–Green tensor:

```math
\gamma_i = \lambda_i^2
```

## Nearly-Incompressible (uncoupled) Ogden model 

```math
\begin{equation}
\label{eq:2}
\Psi = \sum_{i=1}^{N} \frac{c_i}{m_i} (\tilde{\lambda}_1^{m_i}+ \tilde{\lambda}_2^{m_i}+ \tilde{\lambda}_3^{m_i} -3) + \frac{1}{2} K (\log J)^2
\end{equation}
```

!!! note
    **Isochoric (Deviatoric) Stretch Component**

    ```math
    \tilde{\lambda}_i = J^{-1/3} \lambda_i
    ```

Jacobian:

```math
J = \sqrt{\det(\mathbf{C})}
```