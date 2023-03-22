module Gradus

import Base: in
using Base.Threads: @threads
using LinearAlgebra: ×, ⋅, norm, det, dot

using DocStringExtensions
using Parameters

using SciMLBase
using OrdinaryDiffEq
using DiffEqCallbacks
using StaticArrays
using Optim
using DataInterpolations
using VoronoiCells
using FiniteDifferences
using Roots
using ProgressMeter
using Buckets
using QuadGK
using MuladdMacro

using Accessors: @set
using Tullio: @tullio

import ForwardDiff
import GeometryBasics

include("GradusBase/GradusBase.jl")
import .GradusBase:
    E,
    Lz,
    AbstractMetricParameters,
    metric_params,
    metric,
    process_solution,
    process_solution_full,
    GeodesicPoint,
    AbstractGeodesicPoint,
    vector_to_local_sky,
    AbstractMetricParameters,
    geodesic_equation,
    constrain,
    inner_radius,
    metric_type,
    metric_components,
    inverse_metric_components,
    unpack_solution,
    dotproduct,
    propernorm,
    tetradframe,
    lnrbasis,
    lnrframe,
    lowerindices,
    raiseindices,
    StatusCodes,
    AbstractIntegrationParameters,
    IntegrationParameters,
    update_integration_parameters!,
    restrict_ensemble

export AbstractMetricParameters,
    process_solution,
    process_solution_full,
    GeodesicPoint,
    AbstractGeodesicPoint,
    AbstractMetricParameters,
    constrain,
    inner_radius,
    metric_components,
    inverse_metric_components,
    dotproduct,
    propernorm,
    tetradframe,
    lnrbasis,
    lnrframe,
    lowerindices,
    raiseindices,
    StatusCodes,
    AbstractIntegrationParameters,
    IntegrationParameters

# export static arrays too
export SVector, @SVector

"""
    abstract type AbstractPointFunction

Abstract super type for point functions. Must have `f::Function` field.
"""
abstract type AbstractPointFunction end

abstract type AbstractCacheStrategy end
abstract type AbstractRenderCache{M,T} end

abstract type AbstractSkyDomain end
abstract type AbstractGenerator end

"""
    abstract type AbstractAccretionGeometry{T}

Supertype of all accretion geometry. Concrete sub-types must minimally implement
- [`in_nearby_region`](@ref)
- [`has_intersect`](@ref)

Alternativey, for more control, either [`intersects_geometry`](@ref) or [`geometry_collision_callback`](@ref)
may be implemented for a given geometry type.

Geometry intersection calculations are performed by strapping discrete callbacks to the integration
procedure.
"""
abstract type AbstractAccretionGeometry{T} end

"""
    abstract type AbstractAccretionDisc{T} <: AbstractAccretionGeometry{T}

Supertype for axis-symmetric geometry, where certain optimizing assumptions
may be made. Concrete subtypes must implement [`distance_to_disc`](@ref).
"""
abstract type AbstractAccretionDisc{T} <: AbstractAccretionGeometry{T} end

"""
    abstract type AbstractThickAccretionDisc{T} <: AbstractAccretionDisc{T}

Supertype for axis-symmetric geometry that are specified by a height cross-section function. Subtypes are
required to implement [`cross_section`](@ref).
"""
abstract type AbstractThickAccretionDisc{T} <: AbstractAccretionDisc{T} end

"""
    AbstractDiscProfile

Abstract type for binning structures over discs (e.g., radial bins, voronoi).
"""
abstract type AbstractDiscProfile end

abstract type AbstractCoronaModel{T} end

abstract type AbstractDirectionSampler{SkyDomain,Generator} end

abstract type AbstractTraceParameters end

struct EnsembleEndpointThreads end

include("utils.jl")

include("tracing/configuration.jl")
include("tracing/tracing.jl")
include("tracing/geodesic-problem.jl")
include("tracing/constraints.jl")
include("tracing/charts.jl")
include("tracing/callbacks.jl")
include("tracing/utility.jl")
include("tracing/precision-solvers.jl")
include("tracing/radiative-transfer-problem.jl")

include("tracing/method-implementations/auto-diff.jl")

include("image-planes/grids.jl")
include("image-planes/planes.jl")

include("rendering/cache.jl")
include("rendering/rendering.jl")
include("rendering/utility.jl")

include("tracing/method-implementations/first-order.jl")

include("point-functions.jl")

include("orbits/circular-orbits.jl")
include("orbits/orbit-discovery.jl")
include("orbits/orbit-interpolations.jl")

include("accretion-geometry/geometry.jl")
include("accretion-geometry/intersections.jl")
include("accretion-geometry/discs.jl")
include("accretion-geometry/meshes.jl")
include("accretion-geometry/bootstrap.jl")

include("transfer-functions/types.jl")
include("transfer-functions/cunningham-transfer-functions.jl")
include("transfer-functions/integration.jl")

include("corona-to-disc/sky-geometry.jl")
include("corona-to-disc/corona-models.jl")
include("corona-to-disc/disc-profiles.jl")
# needs the types from disc profiles so defer include
include("transfer-functions/transfer-functions-2d.jl")
include("corona-to-disc/flux-calculations.jl")

include("metrics/boyer-lindquist-ad.jl")
include("metrics/boyer-lindquist-fo.jl")
include("metrics/johannsen-ad.jl")
include("metrics/johannsen-psaltis-ad.jl")
include("metrics/morris-thorne-ad.jl")
include("metrics/kerr-refractive-ad.jl")
include("metrics/dilaton-axion-ad.jl")
include("metrics/bumblebee-ad.jl")
include("metrics/kerr-newman-ad.jl")

include("special-radii.jl")
include("redshift.jl")
include("const-point-functions.jl")

include("line-profiles.jl")
include("accretion-geometry/polish-doughnut.jl")

export AbstractPointFunction,
    AbstractCacheStrategy,
    AbstractRenderCache,
    AbstractSkyDomain,
    AbstractGenerator,
    AbstractAccretionGeometry,
    AbstractAccretionDisc,
    AbstractDiscProfile,
    AbstractDirectionSampler

# precompilation help
precompile(
    tracegeodesics,
    (KerrMetric{Float64}, SVector{4,Float64}, SVector{4,Float64}, Tuple{Float64,Float64}),
)
precompile(
    rendergeodesics,
    (KerrMetric{Float64}, SVector{4,Float64}, GeometricThinDisc{Float64}, Float64),
)

end # module
