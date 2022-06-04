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
            # "Tracing"
            # "Rendering"
            "Point Functions" => "overview/point-functions.md",
            # "Callbacks"
            "Available metrics" => "overview/metrics.md"
        ],
        # "Reverberation Lags"
        "Internals" => [
            "Geodesic integration" => "overview/geodesic-integration.md",
            "Implementing new metrics" => "internals/custom-metrics.md",
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
