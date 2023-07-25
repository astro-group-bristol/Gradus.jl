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
import Symbolics

include("GradusBase/GradusBase.jl")
import .GradusBase:
    AbstractCoordinates,
    AbstractTrace,
    AbstractStaticAxisSymmetric,
    BoyerLindquist,
    E,
    Lz,
    AbstractMetric,
    metric,
    unpack_solution,
    unpack_solution_full,
    GeodesicPoint,
    AbstractGeodesicPoint,
    vector_to_local_sky,
    AbstractMetric,
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
    restrict_ensemble,
    _fast_dot

export AbstractMetric,
    AbstractTrace,
    AbstractStaticAxisSymmetric,
    BoyerLindquist,
    AbstractCoordinates,
    unpack_solution,
    unpack_solution_full,
    GeodesicPoint,
    AbstractGeodesicPoint,
    AbstractMetric,
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

abstract type AbstractOpticalProperty end
struct OpticallyThin <: AbstractOpticalProperty end
struct OpticallyThick <: AbstractOpticalProperty end

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
optical_property(::Type{<:AbstractAccretionGeometry}) = OpticallyThick()
optical_property(::T) where {T<:AbstractAccretionGeometry} = optical_property(T)

is_optically_thin(g::AbstractAccretionGeometry) = optical_property(g) isa OpticallyThin

is_finite_disc(::Type{<:AbstractAccretionGeometry}) = true
is_finite_disc(::T) where {T<:AbstractAccretionGeometry} = is_finite_disc(T)

Base.length(::AbstractAccretionGeometry) = 1
Base.iterate(g::AbstractAccretionGeometry) = (g, nothing)
Base.iterate(::AbstractAccretionGeometry, ::Nothing) = nothing

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

struct CompositeGeometry{T,G} <: AbstractAccretionGeometry{T}
    geometry::G
end

"""
    AbstractDiscProfile

Abstract type for binning structures over discs (e.g., radial bins, voronoi).
"""
abstract type AbstractDiscProfile end

"""
    AbstractCoronaModel{T,N}

Type `T`.
"""
abstract type AbstractCoronaModel{T} end

Base.length(::AbstractCoronaModel) = 1
Base.iterate(m::AbstractCoronaModel) = (m, nothing)
Base.iterate(::AbstractCoronaModel, ::Nothing) = nothing

abstract type AbstractDirectionSampler{SkyDomain,Generator} end

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

include("geometry/geometry.jl")
include("geometry/intersections.jl")
include("geometry/discs.jl")
include("geometry/meshes.jl")
include("geometry/composite.jl")
include("geometry/bootstrap.jl")

include("transfer-functions/types.jl")
include("transfer-functions/cunningham-transfer-functions.jl")
include("transfer-functions/integration.jl")

include("corona/samplers.jl")
include("corona/corona-models.jl")
include("corona/disc-profiles.jl")
# needs the types from disc profiles so defer include
include("transfer-functions/transfer-functions-2d.jl")
include("corona/flux-calculations.jl")
include("corona/emissivity.jl")

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
include("geometry/polish-doughnut.jl")

include("plotting-recipes.jl")

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
