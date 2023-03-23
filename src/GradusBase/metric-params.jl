"""
    abstract type AbstractMetric{T} end

Abstract type used to dispatch different geodesic problems.
"""
abstract type AbstractMetric{T} end

Base.length(::AbstractMetric) = 1
Base.iterate(m::AbstractMetric) = (m, nothing)
Base.iterate(::AbstractMetric, ::Nothing) = nothing

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
    constrain(m::AbstractMetric, x, v; μ=0)

Calculate the time component ``v^t`` of a velocity vector `v`, which would constrain the vector at a position `x` as a 
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
constrain(m::AbstractMetric{T}, x, v; μ::T = 0.0) where {T} =
    error("Not implemented for metric parameters $(typeof(m))")

"""
    inner_radius(m::AbstractMetric{T})

Return the innermost valid coordinate relative to the origin, for use in geodesic tracing.

This usually represents some property of the metric, e.g. event horizon radius in Kerr/Schwarzschild metrics, or
throat diameter in worm hole metrics.
"""
inner_radius(::AbstractMetric{T}) where {T} = convert(T, 0.0)

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
