# FerriteHyperelastic

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Aminofa70.github.io/FerriteHyperelastic.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Aminofa70.github.io/FerriteHyperelastic.jl/dev/)
[![Build Status](https://github.com/Aminofa70/FerriteHyperelastic.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Aminofa70/FerriteHyperelastic.jl/actions/workflows/CI.yml?query=branch%3Amain)

# About ComodoFerrite.jl

A Julia Package for finite element analysis and curve fitting of hyperelastic materials.

Built on top of [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/)

This project uses the following Julia packages:

- [Comodo.jl](https://github.com/COMODO-research/Comodo.jl) – for general mesh generation.
- [ComodoFerrite.jl](https://github.com/COMODO-research/ComodoFerrite.jl) – for ferrite simulations and specialized meshes.

Please see the examples in [examples](https://github.com/Aminofa70/FerriteHyperelastic.jl/tree/main/examples)

<img src="https://github.com/Aminofa70/FerriteHyperelastic.jl/blob/main/assets/anim/torsion.gif" alt="ComodoFerrite torsion loading of a cube" width="50%"/>

Supports:


- [x] Plane Stress

- [x] Plane Strain

- [x] 3D

- [x] Curve fitting

To install

```
julia> ]

(@v1.11) pkg> add https://github.com/Aminofa70/FerriteHyperelastic.jl
````





