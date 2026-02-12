# using FerriteHyperelastic
# using Documenter

# DocMeta.setdocmeta!(FerriteHyperelastic, :DocTestSetup, :(using FerriteHyperelastic); recursive=true)

# makedocs(;
#     modules=[FerriteHyperelastic],
#     authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
#     sitename="FerriteHyperelastic.jl",
#     format=Documenter.HTML(;
#         canonical="https://Aminofa70.github.io/FerriteHyperelastic.jl",
#         edit_link="main",
#         assets=String[],
#     ),
#     pages=[
#         "Home" => "index.md",
#     ],
# )

# deploydocs(;
#     repo="github.com/Aminofa70/FerriteHyperelastic.jl",
#     devbranch="main",
# )

using Documenter
using Documenter.Remotes
using FerriteHyperelastic

DocMeta.setdocmeta!(FerriteHyperelastic, :DocTestSetup, :(using FerriteHyperelastic); recursive=true)

makedocs(
    modules=[FerriteHyperelastic],
    sitename="FerriteHyperelastic.jl",
    authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    repo=Remotes.GitHub("Aminofa70", "FerriteHyperelastic.jl"),
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        repolink="https://github.com/Aminofa70/FerriteHyperelastic.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Install" => "install.md",
        "Functions" => "functions.md",
        "Tutorials" => [
            "tutorial curve fitting" => "tutorials/tutorial_curve_fitting.md",
        ]
    ],
)

deploydocs(
    repo = "github.com/Aminofa70/FerriteHyperelastic.jl",
    devbranch = "main",
    push_preview = false,
    versions = [
        "stable" => "v1.0.0",  
        "dev" => "main",
    ],
)

# deploydocs(
#     repo = "github.com/Aminofa70/FerriteHyperelastic.jl",
#     devbranch = "main",
#     push_preview = false,
# )


