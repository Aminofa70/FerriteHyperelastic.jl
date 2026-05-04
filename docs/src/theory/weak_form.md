
## Principle of Minimum Potential Energy

### Derivation of the finite element formulation based on the principle of minimum potential energy.

Total Potential Energy:
```math
\begin{equation}
\label{eq:1}
\mathbf{\Pi} (\mathbf{u}) = U_{int} - U_{ext}
\end{equation}
```

Internal Energy:
```math
\begin{equation}
\label{eq:2}
U_{int} = \int_\Omega \Psi \, d\Omega
\end{equation}
```

External Energy:

::: details

External energy has two contributions: body forces and traction.

There are two types of traction loads:
- **Dead load**: independent of the deformation.
- **Follower load**: depends on the deformation.

:::

:::tabs

== Dead Load

```math
\begin{equation}
\label{eq:3}
    U_{ext} = \oint_{\Gamma_0} \mathbf{t}_0 \mathbf{u}\,d\Gamma_0
\end{equation}
```

== Follower Load 

```math
\begin{equation}
\label{eq:3-3}
U_{ext} = \oint_{\Gamma_0} \mathbf{t} \cdot \mathbf{u}\; J(\mathbf{u}) \sqrt{\mathbf{N} \cdot \mathbf{C}(\mathbf{u})^{-1} \mathbf{N}} \, d\Gamma_0
\end{equation}
```
:::

**Stationarity Condition** at equilibrium, the total potential energy is stationary:

```math
\begin{equation}
\begin{aligned}
\delta \Pi &= 0 \\
\delta U_{\text{int}} - \delta U_{\text{ext}} &= 0 \\
\delta U_{\text{int}} &= \delta U_{\text{ext}}
\end{aligned}
\end{equation}
```
 
The variational internal energy is:

```math
\begin{equation}
\begin{aligned}
\delta U_{\text{int}} &= \int_\Omega \delta \Psi \, d\Omega \\
                      &= \int_\Omega \mathbf{P} : \delta \mathbf{F} \, d\Omega \\
                      &= \int_\Omega \mathbf{P} : \nabla \delta \mathbf{u} \, d\Omega
\end{aligned}
\end{equation}
```

::: info

```math
\mathbf{P} = \frac{\partial \Psi}{\partial \mathbf{F}}
```

```math
\mathbf{S} = 2\,\frac{\partial \Psi}{\partial \mathbf{C}}
```

```math
\mathbf{P} = \mathbf{F}\mathbf{S} \quad \text{(for } \Psi=\Psi(\mathbf{C})\text{)}
```
:::


The variational external energy is:

:::tabs

== Dead Load

```math
\begin{equation}
\label{eq:6}
    \delta U_{ext} = \oint_{\Gamma_0} \mathbf{t}_0 \delta \mathbf{u}\,d\Gamma_0
\end{equation}
```

== Follower Load 
```math
\begin{equation}
\label{eq:6-6}
\delta U_{ext} = \oint_{\Gamma_0} \mathbf{t} \cdot \delta \mathbf{u}\; J(\mathbf{u}) \sqrt{\mathbf{N} \cdot \mathbf{C}(\mathbf{u})^{-1} \mathbf{N}} \, d\Gamma_0
\end{equation}
```

:::