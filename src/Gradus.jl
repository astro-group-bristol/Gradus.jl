module Gradus

import Base: in
using Base.Threads: @threads
using LinearAlgebra: ×, ⋅, norm, det, dot, inv

using Parameters

import DiffEqBase

using SciMLBase
using OrdinaryDiffEq
using DiffEqCallbacks
using StaticArrays
using Optim
using DataInterpolations
using VoronoiCells
using Roots
using ProgressMeter
using Buckets
using QuadGK
using MuladdMacro

using Tullio: @tullio

import ForwardDiff
# todo: temporary fix for https://github.com/JuliaGPU/Metal.jl/issues/229
function Base.sin(x::ForwardDiff.Dual{<:ForwardDiff.Tag{F,Float32}}) where {F}
    s = sin(ForwardDiff.value(x))
    c = cos(ForwardDiff.value(x))
    return ForwardDiff.Dual{typeof(x).parameters[1]}(s, c * ForwardDiff.partials(x))
end
function Base.cos(x::ForwardDiff.Dual{<:ForwardDiff.Tag{F,Float32}}) where {F}
    s = sin(ForwardDiff.value(x))
    c = cos(ForwardDiff.value(x))
    return ForwardDiff.Dual{typeof(x).parameters[1]}(c, -s * ForwardDiff.partials(x))
end
import GeometryBasics
import Symbolics

import FFTW

using EnumX

"""
    StatusCodes.T

Status codes that represent the fate of different geodesics.

- `OutOfDomain`: left the integration chart by travelling to effective infinity (see [`chart_for_metric`](@ref)).
- `WithinInnerBoundary`: left the integration chart by falling into the inner horizon [`inner_radius`](@ref).
- `IntersectedWithGeometry`: geodesics intersected with [`AbstractAccretionGeometry`](@ref) and terminated there. Note
that this status code only applies to optically thick geometry (see [`optical_property`](@ref) for more).
- `NoStatus`: the default status code, which is used to imply no calculation has yet been performed.
"""
@enumx StatusCodes begin
    OutOfDomain
    WithinInnerBoundary
    IntersectedWithGeometry
    NoStatus
end

abstract type AbstractCoordinates end
struct BoyerLindquist{C} <: AbstractCoordinates end
# the following are placeholders
struct KerrSchild{C} <: AbstractCoordinates end
struct FlatCartesian{C} <: AbstractCoordinates end
struct FlatSphericalPolar{C} <: AbstractCoordinates end
struct FlatCylindricalPolar{C} <: AbstractCoordinates end

"""
    abstract type AbstractMetric{T} end

Abstract type used to dispatch different geodesic problems.
"""
abstract type AbstractMetric{T,C} end

Base.length(::AbstractMetric) = 1
Base.iterate(m::AbstractMetric) = (m, nothing)
Base.iterate(::AbstractMetric, ::Nothing) = nothing

# some utility abstract types for dispatching special methods
abstract type AbstractStaticAxisSymmetric{T} <: AbstractMetric{T,BoyerLindquist{(:r, :θ)}} end
abstract type AbstractStaticSphericallySymmetric{T} <: AbstractStaticAxisSymmetric{T} end

"""
    metric_components(m::AbstractMetric{T}, x)

Return a tuple with each non-zero metric component for the metric described by `m` at position
`x`. Note that the position need not be a four-vector, and for specific implementations may
only be a subset of the total manifold coordinates. See specific implementations for subtypes of
[`AbstractMetric`](@ref) for details.
"""
metric_components(m::AbstractMetric, x) = error("Not implemented for metric $(typeof(m))")
inverse_metric_components(m::AbstractMetric, x) =
    error("Not implemented for metric $(typeof(m))")

"""
    geodesic_equation(m::AbstractMetric, x, v)

Calculate the four-acceleration of the geodesic equation for a spacetime given by the metric `m`,
four-position `x` and four-velocity `v`.

A geodesic is the shortest path connecting two points in space. For flat space, this is just a straight line. In
curved space, geodesics are analogous to straight lines between points (e.g. the great circle on a sphere).

The geodesic equation calculates the acceleration experienced by a particle at position ``x^\\mu = (t, r, \\theta, \\phi)`` travelling
with tangential velocity ``v^\\nu = \\text{d} x / \\text{d} \\lambda`` due to the curvature of spacetime. The curvature is calculated from the metric, encoded in the 
[Christoffel symbols](https://en.wikipedia.org/wiki/Christoffel_symbols). The acceleration is then calculated via

```math
\\frac{\\text{d}^2 x^\\mu}{\\text{d} \\lambda^2}
    = - \\Gamma^{\\mu}_{\\phantom{\\mu}\\nu\\sigma}
    \\frac{\\text{d}x^\\nu}{\\text{d} \\lambda}
    \\frac{\\text{d}x^\\sigma}{\\text{d} \\lambda}
```

where ``\\Gamma^{\\mu}_{\\phantom{\\mu}\\nu\\sigma}`` are the Christoffel symbols (of the second kind), and ``\\lambda`` is an affine parameter
that parameterizes the solution.
"""
geodesic_equation(m::AbstractMetric, x, v) =
    error("Not implemented for metric parameters $(typeof(m))")

"""
    _constrain(m::AbstractMetric, x, v; μ=0)

Calculate the time component ``v^t`` of a velocity vector `v`, which would _constrain the vector at a position `x` as a 
geodesic with invariant mass `μ`.

The velocity vector needs to only specify the ``v^r``, ``v^\\theta``, and ``v^\\phi`` component, as the ``v^t`` is constrained in GR by

```math
g_{\\sigma\\nu} v^\\sigma v^\\nu = -\\mu^2,
```

where ``\\mu^2`` is the invariant mass of the particle. This furthermore permits a choice of geodesic to trace. The choices correspond to

- `μ = 0.0` (default): null geodesic
- `μ > 0.0`: time-like geodesic
- `μ < 0.0`: space-like geodesic
"""
_constrain(m::AbstractMetric{T}, x, v; μ::T = 0.0) where {T} =
    error("Not implemented for metric parameters $(typeof(m))")

"""
    inner_radius(m::AbstractMetric{T})

Return the innermost valid coordinate relative to the origin, for use in geodesic tracing.

This usually represents some property of the metric, e.g. event horizon radius in Kerr/Schwarzschild metrics, or
throat diameter in worm hole metrics.
"""
inner_radius(::AbstractMetric{T}) where {T} = zero(T)

"""
    metric_type(m::AbstractMetric{T})

Return the [`AbstractMetric`](@ref) type associated with the metric parameters `m`.
"""
metric_type(m::AbstractMetric) = error("Not implemented for metric parameters $(typeof(m))")


"""
    metric(m::AbstractMetric{T}, u)

Numerically evaluate the corresponding metric for [`AbstractMetric`](@ref), given parameter values `m`
and some point `u`.
"""
metric(m::AbstractMetric, u) = error("Not implemented for metric $(typeof(m))")

restrict_ensemble(::AbstractMetric, ensemble) = ensemble

"""
    AbstractTrace

Parameters that are constant throughout the integration (e.g. mass or frequency) for any
number of geodesics. Also used to dispatch different tracing problems.
"""
abstract type AbstractTrace end

"""
    AbstractIntegrationParameters{M}

Parameters that are made available at each step of the integration, that need not be constant.
For example, the turning points or withing-geometry flags.

The integration parameters should track which spacetime `M` they are parameters for.
Integration parameters must implement
- [`set_status_code!`](@ref)
- [`get_status_code`](@ref)
- [`get_metric`](@ref)

For more complex parameters, may also optionally implement
- [`update_integration_parameters!`](@ref)

See the documentation of each of the above functions for details of their operation.
"""
abstract type AbstractIntegrationParameters{M<:AbstractMetric} end

# type alias, since this is often used
const MutStatusCode = MVector{1,StatusCodes.T}

# TODO: temporary fix for https://github.com/SciML/DiffEqBase.jl/issues/918
function DiffEqBase.anyeltypedual(
    ::AbstractIntegrationParameters{<:AbstractMetric{T}},
) where {T}
    if T <: ForwardDiff.Dual
        T
    else
        Any
    end
end

"""
    update_integration_parameters!(old::AbstractIntegrationParameters, new::AbstractIntegrationParameters)

Update (mutate) the `old` integration parameters to take the value of the `new`. Function should return
the `old`.

Note this function is practically only used to update any mutable fields in the integration parameters,
such as resetting any changes to an original state.
"""
function update_integration_parameters!(
    old::AbstractIntegrationParameters,
    new::AbstractIntegrationParameters,
)
    set_status_code!(old, get_status_code(new))
    old
end

"""
    set_status_code!(p::AbstractIntegrationParameters, status::StatusCodes.T)

Update the status [`StatusCodes`](@ref) in `p` with `status`.
"""
set_status_code!(params::AbstractIntegrationParameters, ::StatusCodes.T) =
    error("Not implemented for $(typeof(params))")

"""
    get_status_code(p::AbstractIntegrationParameters)::StatusCodes.T

Return the status [`StatusCodes`](@ref) in `status`.
"""
get_status_code(params::AbstractIntegrationParameters) =
    error("Not implemented for $(typeof(params))")

"""
    get_metric(p::AbstractIntegrationParameters{M})::M where {M}

Return the [`AbstractMetric`](@ref) `m::M` for which the integration parameters
have been specialised.
"""
get_metric(params::AbstractIntegrationParameters) =
    error("Not implemented for $(typeof(params))")

"""
    abstract type AbstractGeodesicPoint{T}

Supertype for geodesic points, used to store information about specific points along geodesic
trajectories.

!!! note
    Currently limited to storing the start and endpoint of any given trajectory. To keep the
    full geodesic path, it is encouraged to use the `SciMLBase.AbstractODESolution` directly.

Must minimally have the same fields as [`GeodesicPoint`](@ref).
Examples include [`Gradus.FirstOrderGeodesicPoint`](@ref).
"""
abstract type AbstractGeodesicPoint{T} end

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
    abstract type AbstractCoronaModel{T}

The supertype of coronal models, which concrete models must subtype.
Struct implementing `AbstractCoronaModel` must implement minimally
[`sample_position_velocity`](@ref).

For example, adding the lamp-post coronal model
```julia
struct LampPostModel{T} <: AbstractCoronaModel{T}
    height::T
end

function Gradus.sample_position_velocity(m::AbstractMetric, model::LampPostModel)
    # avoid coordinate singularity with a small θ
    x = SVector(0, model.height, 1e-3, 0)
    # ensure velocity is normalized
    g = metric_components(m, SVector(x[2], x[3]))
    v = inv(√(-g[1])) * SVector(1, 0, 0, 0)
    x, v
end
```

Note that [`sample_position_velocity`](@ref) has a number of its own requirements (see that
function's documentation). This function must be implemented as a fallback for other methods.

If special symmetries exist, these may be used in the implementations of higher-order functions, such as
[`emissivity_profile`](@ref). 
"""
abstract type AbstractCoronaModel{T} end

is_point_source(::Type{<:AbstractCoronaModel}) = false
is_point_source(::T) where {T<:AbstractCoronaModel} = is_point_source(T)

Base.length(::AbstractCoronaModel) = 1
Base.iterate(m::AbstractCoronaModel) = (m, nothing)
Base.iterate(::AbstractCoronaModel, ::Nothing) = nothing

"""
abstract type AbstractCoronalSpectrum end

Abstract type specifying interface for the spectra.
"""
abstract type AbstractCoronalSpectrum end

# for broadcasting
Base.length(::AbstractCoronalSpectrum) = 1
Base.iterate(g::AbstractCoronalSpectrum) = (g, nothing)
Base.iterate(::AbstractCoronalSpectrum, ::Nothing) = nothing

coronal_spectrum(spectrum::AbstractCoronalSpectrum, g) =
    error("not implemented for $(typeof(spectrum))")

abstract type AbstractDirectionSampler{SkyDomain,Generator} end

struct EnsembleEndpointThreads end


abstract type AbstractComputationalMethod end
struct TransferFunctionMethod <: AbstractComputationalMethod end
struct BinningMethod <: AbstractComputationalMethod end

include("orthonormalization.jl")
include("utils.jl")
include("solution-processing.jl")

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
include("orbits/orbit-solving.jl")

include("metrics/kerr-metric.jl")
include("metrics/kerr-metric-first-order.jl")
include("metrics/johannsen-ad.jl")
include("metrics/johannsen-psaltis-ad.jl")
include("metrics/morris-thorne-ad.jl")
include("metrics/kerr-refractive-ad.jl")
include("metrics/noz-metric.jl")
include("metrics/dilaton-axion-ad.jl")
include("metrics/bumblebee-ad.jl")
include("metrics/kerr-newman-ad.jl")
include("metrics/kerr-dark-matter.jl")
include("metrics/minkowski.jl")

include("geometry/geometry.jl")
include("geometry/intersections.jl")
include("geometry/discs.jl")
include("geometry/meshes.jl")
include("geometry/composite.jl")
include("geometry/bootstrap.jl")

include("transfer-functions/types.jl")
include("transfer-functions/utils.jl")
include("transfer-functions/cunningham-transfer-functions.jl")
include("transfer-functions/integration.jl")

include("corona/samplers.jl")
include("corona/corona-models.jl")
include("corona/disc-profiles.jl")
# needs the types from disc profiles so defer include
include("transfer-functions/transfer-functions-2d.jl")
include("corona/flux-calculations.jl")
include("corona/emissivity.jl")
include("corona/spectra.jl")

include("special-radii.jl")
include("redshift.jl")
include("const-point-functions.jl")

include("line-profiles.jl")
include("reverberation.jl")

include("plotting-recipes.jl")

if Base.VERSION >= v"1.4.2"
    include("precompile.jl")
end

export AbstractMetric,
    AbstractPointFunction,
    AbstractCacheStrategy,
    AbstractRenderCache,
    AbstractSkyDomain,
    AbstractGenerator,
    AbstractAccretionGeometry,
    AbstractAccretionDisc,
    AbstractDiscProfile,
    AbstractDirectionSampler

export unpack_solution, unpack_solution_full

export metric_components

export StatusCodes

export BinningMethod, TransferFunctionMethod

# export static arrays too
export SVector, @SVector

end # module
