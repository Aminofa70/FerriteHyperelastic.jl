---
title: 'FerriteHyperelastic.jl: A Julia Package for Finite Element Analysis of Hyperelastic Materials'
tags:
  - Julia 
  - Computational Mechanics
  - Finite Element
  - Hyperelasticity 
  - Curve Fitting
authors:
  - name: Amin Alibakhshi
    orcid: 0000-0003-3010-3548
    affiliation: 1  
  - name: Kian Aghani
    orcid: 0000-0002-0707-1157
    affiliation:  2  
  - name: Luis Saucedo-Mora
    orcid: 0000-0003-1091-6875
    affiliation: 1
  - name: Kevin Mattheus Moerman
    orcid: 0000-0003-3768-4269
    affiliation: 3, 
affiliations:
  - name: E.T.S. de Ingeniería Aeronáutica y del Espacio, Universidad Politécnica de Madrid, Pza. Cardenal Cisneros 3, 28040, Madrid, Spain.
    index: 1
  - name: Research assistant, Sahand University of Technology, Tabriz, Iran
    index: 2
  - name: University of Galway, Galway, Ireland
    index: 3
date: 12 Feb 2026
bibliography: paper.bib   
---

# Summary
Finite element (FE) analysis of hyperelastic materials presents unique challenges, particularly in achieving efficient implementation and providing tools for material parameter identification. This article provides an open-source package for FE analysis of hyperelastic materials, built upon the robust [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) [@carlsson2025ferrite] framework in Julia programming language. It supports both two- and three-dimensional FE formulations. A key feature of the package is its integrated curve-fitting module, which employs numerical procedures to determine material constants directly from different experimental data, while also performing stability tests.  Visualization of the geometry, mesh, boundary conditions, and deformation is also provided in the package.

# Statement of need
The efficient implementation of the finite element (FE) for hyperelastic structures presents theoretical and computational challenges [@ahmadi2025plane]. The accurate modeling of hyperelastic materials is vital for engineering applications, making efficient and accessible computational tools crucial for research purposes and industrial applications. The FE software for nonlinear analysis is dominated by two categories: commercial packages and open-source implementations. Commercial software offers robust hyperelasticity modules; however, these solutions present significant barriers, including closed-source architectures and substantial licensing costs. On the other hand, open-source packages are for research purposes and are commonly written in C, C++, or FORTRAN, exemplified by libraries like [FEBio](https://github.com/febiosoftware) [@maas2012febio]. Although these packages offer efficiency and accessibility, they suffer from steep learning curves and verbose syntax that obscures the elegance of the FE procedure.

Computational mechanics software development has experienced advancements with the emergence of high-level computing languages such as Python. Software packages such as [FEniCS](https://fenicsproject.org/) [@baratta2023dolfinx] cover both linear and nonlinear finite element analysis (FEA). However, Python's interpreted nature introduces performance limitations that become an issue in simulations having a high number of degrees of freedom (DOF). This matter becomes more pronounced in element-level operations demanding nested loops.

By the introduction of the Julia programming language in 2012 [@bezanson2017julia], a compelling solution to the "two-language problem" in scientific computing is given by an innovative design combining multiple dispatch, just-in-time compilation, meta-programming capabilities, and optional type annotations. Moreover, the language's syntax closely mirrors mathematical notation, paving the way for direct translation of FE formulations. Within the Julia FE ecosystem, two packages have emerged as foundational frameworks:  [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) [@carlsson2025ferrite] and [Gridap.jl](https://github.com/gridap/Gridap.jl) [@Badia2020].
While both packages provide significant advancements for FEA in Julia, their native support for the domain of hyperelasticity remains limited, which reveals several critical gaps for a researcher or engineer.  Critical gaps include absence of specialized plane stress algorithms, lack of integrated parameter identification tools, and minimal support for common experimental configurations such as uniaxial, biaxial, and pure shear tests.

We present `FerriteHyperelastic`, a specialized Julia package that extends [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) [@carlsson2025ferrite] with comprehensive capabilities for FEA of hyperelastic materials. Our package addresses the aforementioned limitations through several key contributions: (1) integrated curve-fitting modules for material parameter identification and stability check from uniaxial, biaxial, and pure shear experimental data for Neo-Hookean, Mooney-Rivlin, Yeoh, and Ogden models; (2) robust plane stress algorithms employing local Newton iterations with consistent linearization; (3) specialized element formulations (hybrid method) for nearly incompressible materials including mixed methods; (4) visualization tools for pre- and post-processing. The package provides an advancement in accessible, efficient, and mathematically transparent tools for hyperelastic finite element analysis.

# Software design

`FerriteHyperelastic` can be used for the FEA of two- and three-dimensional hyperelastic structures. Thanks to the capabilities of  [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) [@carlsson2025ferrite], in `FerriteHyperelastic`, different mesh types and resolutions can be implemented, Neumann and Dirichlet boundary conditions are specified easily, and results are saved in a VTU format for post-processing. Furthermore, because [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) [@carlsson2025ferrite] supports automatic differentiation, many hyperelastic models can be easily implemented in a small number of lines of code, and related stresses and their differentiation are obtained. It also encompasses the curve-fitting for defining the material parameters in hyperelastic models used in the FEM. 

The package is written in the Julia language and contains different folders. The main folder is src, which includes the module and source functions. The main part of the package is the example folder that supplys different demos for compressible, nearly-incompressible, and fully incompressible problems. To execute the code and perform FEA, the Julia language and `FerriteHyperelastic` should be installed. The package is built on top of other Julia packages, mainly [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) [@carlsson2025ferrite] and others such as [Makie.jl] (https://docs.makie.org/stable/) [@danisch2021makie] (using its backend GLMakie), [GeometryBasics.jl](https://github.com/JuliaGeometry/GeometryBasics.jl) [@Simon_GeometryBasics_2025], and [Tensors.jl](https://github.com/Ferrite-FEM/Tensors.jl) [@carlsson2019tensors]. 


# Acknowledgements

# AI usage disclosure


# References