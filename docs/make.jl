using ShrodingerSolver
using Documenter

DocMeta.setdocmeta!(ShrodingerSolver, :DocTestSetup, :(using ShrodingerSolver); recursive=true)

makedocs(;
    modules=[ShrodingerSolver],
    authors="walexaindre <walexaindre@hotmail.com
> and contributors",
    repo="https://github.com/walexaindre/ShrodingerSolver.jl/blob/{commit}{path}#{line}",
    sitename="ShrodingerSolver.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://walexaindre.github.io/ShrodingerSolver.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/walexaindre/ShrodingerSolver.jl",
    devbranch="main",
)
