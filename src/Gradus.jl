module Gradus

include("GradusBase/GradusBase.jl")
include("GeodesicTracer/GeodesicTracer.jl")
include("FirstOrderMethods/FirstOrderMethods.jl")
include("Rendering/Rendering.jl")
include("AccretionGeometry/AccretionGeometry.jl")
include("DiscProfiles/DiscProfiles.jl")

using .GradusBase
using .GeodesicTracer
using .FirstOrderMethods
using .Rendering
using .AccretionGeometry
using .DiscProfiles

using StaticArrays
using Parameters
using DocStringExtensions

# GradusBase
export AbstractMetricParams,
    metric_params,
    metric,
    get_endpoint,
    get_startpoint,
    get_point,
    GeodesicPoint,
    vector_to_local_sky,
    AbstractMetricParams,
    geodesic_eq,
    geodesic_eq!,
    constrain,
    on_chart,
    inner_radius,
    metric_type

# GeodesicTracer
export tracegeodesics,
    map_impact_parameters,
    AbstractAutoDiffMetricParams,
    AbstractAutoDiffStaticAxisSymmetricParams,
    metric_components

# FirstOrderMethods
export AbstractFirstOrderMetricParams, FirstOrderGeodesicPoint, BoyerLindquistFO

# Rendering
export rendergeodesics,
    prerendergeodesics, PointFunction, FilterPointFunction, ConstPointFunctions, apply

# AccretionGeometry
export AbstractAccretionGeometry,
    AbstractAccretionDisc, MeshAccretionGeometry, GeometricThinDisc, in_polygon, get_area

# DiscProfiles
export AbstractCoronaModel,
    LampPostModel,
    renderprofile,
    LowerHemisphere,
    BothHemispheres,
    EvenSampler,
    WeierstrassSampler,
    RandomGenerator,
    GoldenSpiralGenerator,
    VoronoiDiscProfile

# pre-defined metrics
include("metrics/metrics.jl")

export BoyerLindquistAD, BoyerLindquistFO, JohannsenAD, MorrisThorneAD

# downstream modules and work
include("special-radii.jl")

include("AccretionFormulae/AccretionFormulae.jl")

using .AccretionFormulae
export solve_equitorial_circular_orbit, trace_equitorial_circular_orbit

end # module
