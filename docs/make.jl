using ParametricSurfaces
using Documenter

DocMeta.setdocmeta!(ParametricSurfaces, :DocTestSetup, :(using ParametricSurfaces); recursive=true)

makedocs(;
    modules=[ParametricSurfaces],
    authors="Luiz M. Faria <maltezfaria@gmail.com> and contributors",
    repo="https://github.com/WaveProp/ParametricSurfaces.jl/blob/{commit}{path}#{line}",
    sitename="ParametricSurfaces.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://WaveProp.github.io/ParametricSurfaces.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/WaveProp/ParametricSurfaces.jl",
    devbranch="main",
)
