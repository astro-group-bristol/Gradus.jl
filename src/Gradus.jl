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


# GradusBase
export AbstractMetricParams,
    metric_params,
    metric,
    get_endpoint,
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
    prerendergeodesics, PointFunction, FilterPointFunction, ConstPointFunctions

# AccretionGeometry
export AbstractAccretionGeometry,
    AbstractAccretionDisc, MeshAccretionGeometry, GeometricThinDisc

# DiscProfiles
export AbstractCoronaModel,
    LampPostModel,
    tracegeodesics,
    renderprofile,
    LowerHemisphere,
    BothHemispheres,
    RandomSampler,
    EvenSampler,
    WeierstrassSampler


end # module
