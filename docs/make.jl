using FerriteHyperelastic
using Documenter

DocMeta.setdocmeta!(FerriteHyperelastic, :DocTestSetup, :(using FerriteHyperelastic); recursive=true)

makedocs(;
    modules=[FerriteHyperelastic],
    authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    sitename="FerriteHyperelastic.jl",
    format=Documenter.HTML(;
        canonical="https://Aminofa70.github.io/FerriteHyperelastic.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Aminofa70/FerriteHyperelastic.jl",
    devbranch="main",
)
