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

import ForwardDiff
import Roots

# GradusBase
export AbstractMetricParams,
    AbstractGeodesicPoint,
    metric_params,
    metric,
    getgeodesicpoint,
    GeodesicPoint,
    vector_to_local_sky,
    AbstractMetricParams,
    geodesic_eq,
    geodesic_eq!,
    constrain,
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
    prerendergeodesics, PointFunction, FilterPointFunction, apply

# AccretionGeometry
export AbstractAccretionGeometry,
    AbstractAccretionDisc,
    MeshAccretionGeometry,
    GeometricThinDisc,
    inpolygon,
    getarea,
    CircularOrbits

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
    VoronoiDiscProfile,
    findindex,
    getareas,
    getproperarea,
    bin_transfer_function

# pre-defined metrics
include("metrics/metrics.jl")

export BoyerLindquistAD,
    BoyerLindquistFO, JohannsenAD, MorrisThorneAD, KerrRefractiveAD, DilatonAxionAD

# downstream modules and work
include("special-radii.jl")

include("AccretionFormulae/AccretionFormulae.jl")

using .AccretionFormulae
export solve_equitorial_circular_orbit,
    trace_equitorial_circular_orbit,
    interpolate_plunging_velocities,
    PlungingInterpolation,
    interpolate_redshift

include("const-point-functions.jl")
export ConstPointFunctions


end # module
