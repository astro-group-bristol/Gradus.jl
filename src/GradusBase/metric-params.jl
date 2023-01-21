"""
    abstract type AbstractMetricParams{T} end

Abstract type used to dispatch different geodesic problems.
"""
abstract type AbstractMetricParams{T} end


# contains the full metric components (this type needed for DiffGeoSymbolics)
abstract type AbstractMetric{T} <: AbstractMatrix{T} end

metric_params(m::AbstractMetric{T}) where {T} =
    error("Not implemented for metric $(typeof(m))")

Base.length(::AbstractMetricParams) = 1
Base.iterate(m::AbstractMetricParams) = (m, nothing)
Base.iterate(::AbstractMetricParams, ::Nothing) = nothing

"""
    metric_components(m::AbstractMetricParams{T}, x)

Return a tuple with each non-zero metric component for the metric described by `m` at position
`x`. Note that the position need not be a four-vector, and for specific implementations may
only be a subset of the total manifold coordinates. See specific implementations for subtypes of
[`AbstractMetricParams`](@ref) for details.
"""
metric_components(m::AbstractMetricParams, x) =
    error("Not implemented for metric $(typeof(m))")
inverse_metric_components(m::AbstractMetricParams, x) =
    error("Not implemented for metric $(typeof(m))")

"""
    geodesic_eq(m::AbstractMetricParams{T}, u, v)
    geodesic_eq!(m::AbstractMetricParams{T}, u, v)

Calculate the acceleration components of the geodesic equation given a position `u`, a velocity `v`, and a metric `m`.
"""
geodesic_eq(m::AbstractMetricParams, u, v) =
    error("Not implemented for metric parameters $(typeof(m))")
geodesic_eq!(m::AbstractMetricParams, u, v) =
    error("Not implemented for metric parameters $(typeof(m))")

"""
    constrain(m::AbstractMetricParams{T}, u, v; μ::T=0.0)

Calculate time component ``v^t`` which would constrain a velocity vector `v` at position `x`
as a geodesic with mass `μ`.
"""
constrain(m::AbstractMetricParams{T}, u, v; μ::T = 0.0) where {T} =
    error("Not implemented for metric parameters $(typeof(m))")

"""
    inner_radius(m::AbstractMetricParams{T})

Return the innermost valid coordinate relative to the origin, for use in geodesic tracing.

This usually represents some property of the metric, e.g. event horizon radius in Kerr/Schwarzschild metrics, or
throat diameter in worm hole metrics.
"""
inner_radius(::AbstractMetricParams{T}) where {T} = convert(T, 0.0)

"""
    metric_type(m::AbstractMetricParams{T})

Return the [`AbstractMetric`](@ref) type associated with the metric parameters `m`.
"""
metric_type(m::AbstractMetricParams) =
    error("Not implemented for metric parameters $(typeof(m))")


"""
    metric(m::AbstractMetricParams{T}, u)

Numerically evaluate the corresponding metric for [`AbstractMetricParams`](@ref), given parameter values `m`
and some point `u`.
"""
metric(m::AbstractMetricParams, u) = error("Not implemented for metric $(typeof(m))")
