"""
    abstract type AbstractMetricParameters{T} end

Abstract type used to dispatch different geodesic problems.
"""
abstract type AbstractMetricParameters{T} end


# contains the full metric components (this type needed for DiffGeoSymbolics)
abstract type AbstractMetric{T} <: AbstractMatrix{T} end

metric_params(m::AbstractMetric{T}) where {T} =
    error("Not implemented for metric $(typeof(m))")

Base.length(::AbstractMetricParameters) = 1
Base.iterate(m::AbstractMetricParameters) = (m, nothing)
Base.iterate(::AbstractMetricParameters, ::Nothing) = nothing

"""
    metric_components(m::AbstractMetricParameters{T}, x)

Return a tuple with each non-zero metric component for the metric described by `m` at position
`x`. Note that the position need not be a four-vector, and for specific implementations may
only be a subset of the total manifold coordinates. See specific implementations for subtypes of
[`AbstractMetricParameters`](@ref) for details.
"""
metric_components(m::AbstractMetricParameters, x) =
    error("Not implemented for metric $(typeof(m))")
inverse_metric_components(m::AbstractMetricParameters, x) =
    error("Not implemented for metric $(typeof(m))")

raw"""
    geodesic_equation(m::AbstractMetricParameters, x, v)

Calculate the four-acceleration of the geodesic equation for a spacetime given by the metric `m`,
four-position `x` and four-velocity `v`.

A geodesic is the shortest path connecting two points in space. For flat space, this is just a straight line. In
curved space, geodesics are analogous to straight lines between points (e.g. the great circle on a sphere).

The geodesic equation calculates the acceleration experienced by a particle at position ``x^\mu = (t, r, \theta, \phi)`` travelling
with tangential velocity ``v^\nu = \text{d} x / \text{d} \lambda`` due to the curvature of spacetime. The curvature is calculated from the metric, encoded in the 
[Christoffel symbols](https://en.wikipedia.org/wiki/Christoffel_symbols). The acceleration is then calculated via

```math
\frac{\text{d}^2 x^\mu}{\text{d} \lambda^2}
    = - \Gamma^{\mu}_{\phantom{\mu}\nu\sigma}
    \frac{\text{d}x^\nu}{\text{d} \lambda}
    \frac{\text{d}x^\sigma}{\text{d} \lambda}
```

where ``\Gamma^{\mu}_{\phantom{\mu}\nu\sigma}`` are the Christoffel symbols (of the second kind), and ``\lambda`` is an affine parameter
that parameterizes the solution.
"""
geodesic_equation(m::AbstractMetricParameters, x, v) =
    error("Not implemented for metric parameters $(typeof(m))")

"""
    constrain(m::AbstractMetricParameters{T}, u, v; μ::T=0.0)

Calculate time component ``v^t`` which would constrain a velocity vector `v` at position `x`
as a geodesic with mass `μ`.
"""
constrain(m::AbstractMetricParameters{T}, u, v; μ::T = 0.0) where {T} =
    error("Not implemented for metric parameters $(typeof(m))")

"""
    inner_radius(m::AbstractMetricParameters{T})

Return the innermost valid coordinate relative to the origin, for use in geodesic tracing.

This usually represents some property of the metric, e.g. event horizon radius in Kerr/Schwarzschild metrics, or
throat diameter in worm hole metrics.
"""
inner_radius(::AbstractMetricParameters{T}) where {T} = convert(T, 0.0)

"""
    metric_type(m::AbstractMetricParameters{T})

Return the [`AbstractMetric`](@ref) type associated with the metric parameters `m`.
"""
metric_type(m::AbstractMetricParameters) =
    error("Not implemented for metric parameters $(typeof(m))")


"""
    metric(m::AbstractMetricParameters{T}, u)

Numerically evaluate the corresponding metric for [`AbstractMetricParameters`](@ref), given parameter values `m`
and some point `u`.
"""
metric(m::AbstractMetricParameters, u) = error("Not implemented for metric $(typeof(m))")

restrict_ensemble(::AbstractMetricParameters, ensemble) = ensemble
