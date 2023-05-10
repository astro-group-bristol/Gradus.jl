push!(LOAD_PATH, "src")

using Documenter
using Gradus



makedocs(
    modules = [Gradus],
    clean = true,
    sitename = "Gradus.jl Documentation",
    pages = [
        "Home" => "index.md",
        "Overview" => [
            "Getting started" => "getting-started.md",
            "Problems and solvers" => "problems-and-solvers.md",
            "Point functions" => "overview/point-functions.md",
            "Catalogue of metrics" => "overview/metrics.md",
            "Accretion geometry" => "overview/accretion-geometry.md",
            "Implementing new metrics" => "internals/custom-metrics.md",
            "Examples" => "examples/examples.md",
        ],
        "Advanced" => [
            "Geodesic integration" => "overview/geodesic-integration.md",
            "Custom tracing" => "internals/custom-traces.md",
            "Parallelism and ensembles" => "internals/parallelism.md",
            "Special radii" => "internals/special-radii.md",
        ],
        "API" =>
            [
                "Gradus" => "api-documentation/Gradus.md",
                "GradusBase" => "api-documentation/GradusBase.md",
            ] |> sort,
    ],
)

deploydocs(repo = "github.com/astro-group-bristol/Gradus.jl.git")
