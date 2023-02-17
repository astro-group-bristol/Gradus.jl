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
            "Available metrics" => "overview/metrics.md",
            "Accretion geometry" => "overview/accretion-geometry.md",
            "Examples" => "examples/examples.md",
        ],
        "Internals" => [
            "Geodesic integration" => "overview/geodesic-integration.md",
            "Implementing new metrics" => "internals/custom-metrics.md",
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
