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
    modules = [FerriteHyperelastic],
    sitename = "FerriteHyperelastic.jl",
    authors = "Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    repo = Remotes.GitHub("Aminofa70", "FerriteHyperelastic.jl"),
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        repolink = "https://github.com/Aminofa70/FerriteHyperelastic.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Install" => "install.md",
        "Functions" => "functions.md",
        "Instruction" => "instruction.md",
        "Examples" => [
            "plane strain tension" => "examples/plane_strain_tension.md",
        ],
    ],
)

# deploydocs(
#     repo = "github.com/Aminofa70/FerriteHyperelastic.jl",
#     devbranch = "main",
#     push_preview = false,
# )

deploydocs(
    repo = "https://github.com/username/your-repo.git",  # your GitHub repo
    branch = "gh-pages",                                   # branch for docs
    make = true                                            # build docs before deploying
)
