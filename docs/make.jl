using FerriteHyperelastic
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(FerriteHyperelastic, :DocTestSetup, :(using FerriteHyperelastic); recursive=true)

makedocs(
    modules=[FerriteHyperelastic],
    sitename="FerriteHyperelastic.jl",
    authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    format=DocumenterVitepress.MarkdownVitepress(
        repo="github.com/Aminofa70/FerriteHyperelastic.jl",
        devbranch="main",
        devurl="dev",
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => [
            "Hyperelastic Models" => [
                "Neo-Hookean Model" => "theory/neo_hookean.md",
                "Mooney-Rivlin Model" => "theory/mooney_rivlin.md",
                "Ogden Model" => "theory/ogden.md",
            ],
            "FEM Theory" => [
                "Weak Formulation" => "theory/weak_form.md",
               # "Finite Element Discretization" => "theory/fem/discretization.md",
            ],
        ],
        "Curve Fitting" => [
            "Instruction Curve Fitting" => "tutorials/tutorial_curve_fitting.md",
        ],
        "FEM Ferrite Solver" => [
            "2D" => [
                "Uniaxial Plane Stress Displacement" => "2D/demo_0001_uniaxial_plane_stress_displacement.md",
                "Uniaxial Plane Strain Traction" => "2D/demo_0001_uniaxial_plane_stress_displacement.md",
            ],
            "3D" => [
                "Compressible" => [
                ],
                "Nearly-Incompressible" => [
                ],
                "Incompressible" => [
                ],
            ],
        ],
        "FEM FESolvers Package" => [
            "demo_0001_uniaxial_plane_stress_displacement " => "2D/demo_0001_uniaxial_plane_stress_displacement.md"
        ],
        "API Reference" => "api.md",
    ]
)

DocumenterVitepress.deploydocs(;
    repo="github.com/Aminofa70/FerriteHyperelastic.jl.git",
    target=joinpath(@__DIR__, "build"),
    branch="gh-pages",
    devbranch="main",
    push_preview=true
)