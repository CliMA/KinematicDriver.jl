# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

pages = Any[
    "Home" => "index.md",
    "Simulation setups" => "simulation_setups.md",
    "Calibration features" => "calibration_features.md",
    "References" => "References.md",
]

mathengine = MathJax(Dict(:TeX => Dict(:equationNumbers => Dict(:autoNumber => "AMS"), :Macros => Dict())))

format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", mathengine = mathengine, collapselevel = 1)

makedocs(
    bib,
    sitename = "Kinematic1D.jl",
    strict = true,
    format = format,
    checkdocs = :exports,
    clean = true,
    doctest = true,
    #modules = [Kinematic1D],
    pages = pages,
)

deploydocs(
    repo = "github.com/CliMA/Kinematic1D.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
