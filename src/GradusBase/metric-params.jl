"""
    abstract type AbstractMetricParams{T} end

Abstract type used to dispatch different geodesic problems.
"""
abstract type AbstractMetricParams{T} end


# contains the full metric components (this type needed for DiffGeoSymbolics)
abstract type AbstractMetric{T} <: AbstractMatrix{T} end

metric_params(m::AbstractMetric{T}) where {T} =
    error("Not implemented for metric $(typeof(m))")

"""
    geodesic_eq(m::AbstractMetricParams{T}, u, v)
    geodesic_eq!(m::AbstractMetricParams{T}, u, v)

Calculate the acceleration components of the geodesic equation given a position `u`, a velocity `v`, and a metric `m`.
"""
geodesic_eq(m::AbstractMetricParams{T}, u, v) where {T} =
    error("Not implemented for metric parameters $(typeof(m))")
geodesic_eq!(m::AbstractMetricParams{T}, u, v) where {T} =
    error("Not implemented for metric parameters $(typeof(m))")

"""
    constrain(m::AbstractMetricParams{T}, u, v; μ::T=0.0f0)

Give time component which would constrain a velocity vector `v` at position `x` to be a
geodesic with mass `μ`.
"""
constrain(m::AbstractMetricParams{T}, u, v; μ::T = 0.0) where {T} =
    error("Not implemented for metric parameters $(typeof(m))")

"""
    on_chart(m::AbstractMetricParams{T}, u)

Check if point `u` is a valid point for the metric described by `m`.

Returns false is `u` is a singularity.
"""
on_chart(m::AbstractMetricParams{T}, u) where {T} = !(sum(u) ≈ 0)


"""
    inner_radius(m::AbstractMetricParams{T})

Return the innermost valid coordinate relative to the origin, for use in geodesic tracing.

This usually represents some property of the metric, e.g. event horizon radius in Kerr/Schwarzschild metrics, or 
throat diameter in worm hole metrics.
"""
inner_radius(m::AbstractMetricParams{T}) where {T} = convert(T, 0.0)

"""
    metric_type(m::AbstractMetricParams{T})

Return the [`AbstractMetric`](@ref) type associated with the metric parameters `m`.
"""
metric_type(m::AbstractMetricParams{T}) where {T} =
    error("Not implemented for metric parameters $(typeof(m))")


"""
    metric(m::AbstractMetricParams{T}, u) 

Numerically evaluate the corresponding metric for [`AbstractMetricParams`](@ref), given parameter values `m`
and some point `u`.
"""
metric(m::AbstractMetricParams{T}, u) where {T} =
    error("Not implemented for metric $(typeof(m))")

# do we actually want to support this?
# since if it's a symbolic matrix, you can subs other ways better?
#"""
#    metric(m::AbstractMetric{T}, u) 
#
#Evaluate the metric at a point `u`.
#"""
#metric(m::AbstractMetric{T}, u) where {T} = error("Not implemented for metric $(typeof(m))")

export AbstractMetricParams,
    geodesic_eq, geodesic_eq!, constrain, on_chart, inner_radius, metric_type
