module ConstPointFunctions
import ..Rendering: PointFunction, FilterPointFunction
import ..AccretionFormulae: _redshift_guard

const filter_early_term =
    FilterPointFunction((m, gp, max_time; kwargs...) -> gp.t2 < max_time, NaN)

const filter_intersected =
    FilterPointFunction((m, gp, max_time; kwargs...) -> gp.retcode == :Intersected, NaN)

const affine_time = PointFunction((m, gp, max_time; kwargs...) -> gp.t2)

const shadow = affine_time ∘ filter_early_term

"""
    redshift(m::AbstractMetricParams, gp::GeodesicPoint, max_time)

Calculate the analytic redshift at a given geodesic point, assuming equitorial, geometrically
thin accretion disc. Implementation depends on the metric type. Currently implemented for

- [`BoyerLindquistAD`](@ref)
- [`BoyerLindquistFO`](@ref)

# Notes

Wraps calls to [`AccretionFormulae._redshift_guard`](@ref) to dispatch different implementations.
"""
const redshift = PointFunction(_redshift_guard)

end # module
