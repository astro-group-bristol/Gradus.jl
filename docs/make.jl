push!(LOAD_PATH,"src")

using Documenter
using Gradus

makedocs(
    modules=[Gradus],
    clean=true,
    sitename="Gradus.jl Documentation",

    pages = [
        "Home" => "index.md",
        "Overview" => [
            # "Examples"
            "Geodesic integration" => "overview/geodesic-integration.md"
            "Implemented metrics" => "overview/metrics.md"
        ],
        "Internals" => [
            "Custom metrics" => "internals/custom-metrics.md"
        ],
        "Module API" => [
            "Gradus" => "api-documentation/Gradus.md",
            "GradusBase" => "api-documentation/GradusBase.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/astro-group-bristol/Gradus.jl.git"
)
