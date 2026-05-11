
## Principle of Minimum Potential Energy

### Weak Form Formulation

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
\label{eq:7}
    \delta U_{ext} = \oint_{\Gamma_0} \mathbf{t}_0 \delta \mathbf{u}\,d\Gamma_0
\end{equation}
```

== Follower Load 
```math
\begin{equation}
\label{eq:7-7}
\delta U_{ext} = \oint_{\Gamma_0} \mathbf{t} \cdot \delta \mathbf{u}\; J(\mathbf{u}) \sqrt{\mathbf{N} \cdot \mathbf{C}(\mathbf{u})^{-1} \mathbf{N}} \, d\Gamma_0
\end{equation}
```

:::

The residual is defined as:

:::tabs

== Dead Load
```math
R(\mathbf{u}) = \int_\Omega \mathbf{P} : \nabla \delta \mathbf{u} \, d\Omega- \oint_{\Gamma_0} \mathbf{t}_0 \delta \mathbf{u}\,d\Gamma_0
```

== Follower Load 

```math
R(\mathbf{u}) = \int_\Omega \mathbf{P} : \nabla \delta \mathbf{u} \, d\Omega - \oint_{\Gamma_0} \mathbf{t} \cdot \delta \mathbf{u}\; J(\mathbf{u}) \sqrt{\mathbf{N} \cdot \mathbf{C}(\mathbf{u})^{-1} \mathbf{N}} \, d\Gamma_0
```
:::

## Linearization of residual 

::: info

Taylor series expansion is 
```math
f(x) = f(a) + f'(a)(x-a) + \frac{1}{2} f''(a)(x-a)^2 + ...
```
For linearization of residual, we need forward shift form of the Taylor series.
For the general form,

```math
f(x + h) = f(x) + h f'(x) + \frac{h^2}{2!}f''(x) + ...
```
::: 

In the finite element, we have
```math
u^{k+1} = u^{k} + \Delta u 
```

Then, Taylor is (Linearization)
```math
R(u^{k+1}) = R(u^{k} + \Delta u) \approx R(u^{k}) + R'(u^{k}) \Delta u
```
The tangent stiffness matrix 
```math
K(u^k) := R'(u^{k}) = \frac{\partial R}{\partial u^k}
```
Then, final form of the linearization is
```math
R(u^{k+1}) = R(u^{k} + \Delta u) \approx R(u^{k}) + K \Delta u
```
Equilibrium enforces
```math
R(u^{k+1}) = R(u^{k} + \Delta u) = 0
```
Then
```math
K(u^k) \Delta u = - R(u)
```
For the discretization of the linearized equilibrium equations and the residual, please refer to [Finite Element Discretization](discretization.md)


