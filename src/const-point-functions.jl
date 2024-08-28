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
    AbstractMetric,
    interpolate_redshift
# for doc bindings
import ..Gradus

"""
    filter_early_term(m::AbstractMetric, gp::AbstractGeodesicPoint, max_time)

A [`FilterPointFunction`](@ref) that filters geodesics that termined early (i.e., did not reach maximum integration time or effective infinity).
Default: `NaN`.
"""
function filter_early_term(T::Type = Float64)
    FilterPointFunction((m, gp, max_time; kwargs...) -> gp.λ_max < max_time, T(NaN))
end

"""
    filter_intersected(m::AbstractMetric, gp::AbstractGeodesicPoint, max_time)

A [`FilterPointFunction`](@ref) that filters geodesics which intersected with the accretion
disc. Default: `NaN`.
"""
function filter_intersected(T::Type = Float64)
    FilterPointFunction(
        (m, gp, max_time; kwargs...) -> status(gp) == StatusCodes.IntersectedWithGeometry,
        T(NaN),
    )
end

"""
    affine_time(m::AbstractMetric, gp::AbstractGeodesicPoint, max_time)

A [`PointFunction`](@ref) returning the affine integration time at the endpoint of the geodesic.
"""
function affine_time()
    PointFunction((m, gp, max_time; kwargs...) -> gp.λ_max)
end

"""
    shadow(m::AbstractMetric, gp::AbstractGeodesicPoint, max_time)

A [`PointFunction`](@ref) which colours the shadow of the black hole for any disc-less render.
Equivalent to `ConstPointFunctions.affine_time ∘ ConstPointFunctions.filter_early_term`.
"""
function shadow(T::Type = Float64)
    affine_time() ∘ filter_early_term(T)
end

"""
    redshift(m::AbstractMetric)

Returns a [`PointFunction`](@ref).

Calculate the analytic redshift at a given geodesic point, assuming equatorial, geometrically
thin accretion disc. Implementation depends on the metric type. Currently implemented for

- [`Gradus.KerrMetric`](@ref)
- [`Gradus.KerrSpacetimeFirstOrder`](@ref)

# Notes

Wraps calls to [`Gradus._redshift_guard`](@ref) to dispatch different implementations.
"""
redshift(::KerrMetric, _) = PointFunction(_redshift_guard)
redshift(::KerrSpacetimeFirstOrder, _) = PointFunction(_redshift_guard)
redshift(m::AbstractMetric, u; kwargs...) = interpolate_redshift(m, u; kwargs...)

end # module

export ConstPointFunctions
