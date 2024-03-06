push!(LOAD_PATH, "src")

using Documenter
using Gradus



makedocs(
    modules = [Gradus],
    clean = true,
    sitename = "Gradus.jl Documentation",
    warnonly = true,
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
        "Examples" => "examples.md",
        "Reference & walkthroughs" => [
            "Catalogue of metrics" => "metrics.md",
            "Accretion geometry" => "accretion-geometry.md",
            "Problems and solvers" => "problems-and-solvers.md",
            "Point functions" => "point-functions.md",
            "Disc emissivity" => "emissivity.md",
        ],
        "Geodesics and integration" => [
            "Implementing new metrics" => "custom-metrics.md",
            "Geodesic integration" => "geodesic-integration.md",
            "Parallelism and ensembles" => "parallelism.md",
            "Special radii" => "special-radii.md",
        ],
        "API" => ["Gradus" => "api-documentation/Gradus.md"] |> sort,
    ],
)

deploydocs(repo = "github.com/astro-group-bristol/Gradus.jl.git")
