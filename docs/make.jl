push!(LOAD_PATH,joinpath(@__DIR__, ".."))
using Documenter, ParametricSurfaces

makedocs(
    modules = [ParametricSurfaces],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Luiz M. Faria",
    sitename = "ParametricSurfaces.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/maltezfaria/ParametricSurfaces.jl.git",
)
