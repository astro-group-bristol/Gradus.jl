"""
    module ConstPointFunctions

Module defining a number of `const` [`Gradus.AbstractPointFunction`](@ref), serving different utility
or common purposes for analysis.
"""
module ConstPointFunctions
using ..Gradus:
    PointFunction,
    FilterPointFunction,
    _redshift_guard,
    StatusCodes,
    KerrMetric,
    KerrSpacetimeFirstOrder,
    AbstractMetricParams,
    interpolate_redshift
# for doc bindings
import ..Gradus

"""
    filter_early_term(m::AbstractMetricParams, gp::AbstractGeodesicPoint, max_time)

A [`FilterPointFunction`](@ref) that filters geodesics that termined early (i.e., did not reach maximum integration time or effective infinity).
Default: `NaN`.
"""
const filter_early_term =
    FilterPointFunction((m, gp, max_time; kwargs...) -> gp.t2 < max_time, NaN)

"""
    filter_intersected(m::AbstractMetricParams, gp::AbstractGeodesicPoint, max_time)

A [`FilterPointFunction`](@ref) that filters geodesics which intersected with the accretion
disc. Default: `NaN`.
"""
const filter_intersected = FilterPointFunction(
    (m, gp, max_time; kwargs...) -> gp.status == StatusCodes.IntersectedWithGeometry,
    NaN,
)

"""
    affine_time(m::AbstractMetricParams, gp::AbstractGeodesicPoint, max_time)

A [`PointFunction`](@ref) returning the affine integration time at the endpoint of the geodesic.
"""
const affine_time = PointFunction((m, gp, max_time; kwargs...) -> gp.t2)

"""
    shadow(m::AbstractMetricParams, gp::AbstractGeodesicPoint, max_time)

A [`PointFunction`](@ref) which colours the shadow of the black hole for any disc-less render.
Equivalent to `ConstPointFunctions.affine_time ∘ ConstPointFunctions.filter_early_term`.
"""
const shadow = affine_time ∘ filter_early_term

"""
    redshift(m::AbstractMetricParams)

Returns a [`PointFunction`](@ref).

Calculate the analytic redshift at a given geodesic point, assuming equitorial, geometrically
thin accretion disc. Implementation depends on the metric type. Currently implemented for

- [`Gradus.KerrMetric`](@ref)
- [`Gradus.KerrSpacetimeFirstOrder`](@ref)

# Notes

Wraps calls to [`Gradus._redshift_guard`](@ref) to dispatch different implementations.
"""
redshift(::KerrMetric, _) = PointFunction(_redshift_guard)
redshift(::KerrSpacetimeFirstOrder, _) = PointFunction(_redshift_guard)
redshift(m::AbstractMetricParams, u; kwargs...) = interpolate_redshift(m, u; kwargs...)

end # module

export ConstPointFunctions
