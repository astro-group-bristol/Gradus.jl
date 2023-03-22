"""
    $(TYPEDSIGNATURES)

Innermost stable circular orbit (ISCO), defined by
```math
    \\frac{\\text{d}}{\\text{d}r} \\left( \\frac{E}{\\mu} \\right) = 0.
```
Uses analytic solutions if known for that metric, else uses a root finder to calculate
the radius at which the above condition is met.
"""
isco(m::AbstractMetricParameters) = error("Not implemented for $(typeof(m)).")

function isco(m::AbstractStaticAxisSymmetricParameters, lower_bound, upper_bound; kwargs...)
    if lower_bound == upper_bound
        error(
            "No boundaries for minimization could be determined. It is likely this configuration does not have an ISCO solution.",
        )
    end
    dE(r) = ForwardDiff.derivative(x -> CircularOrbits.energy(m, x; kwargs...), r)
    d2E(r) = ForwardDiff.derivative(dE, r)

    Roots.find_zero((dE, d2E), (lower_bound, upper_bound))
end

function isco(
    m::AbstractStaticAxisSymmetricParameters{T};
    max_upper_bound = T(100),
    step = T(0.005),
    kwargs...,
) where {T}
    lower_bound, upper_bound = find_isco_bounds(m, max_upper_bound, step; kwargs...)
    isco(m, lower_bound, upper_bound, ; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Brute force method for finding the first radius at which ``\\frac{E}{\\mu} > 1``. The method
calculates ``E`` with [`CircularOrbits.energy`](@ref), where `r` steps from `upper_bound` to
``r = r_\\text{g}`` in steps of `step`.

Returns `T(0.0)` if no such radius found.
"""
function find_isco_bounds(
    m::AbstractStaticAxisSymmetricParameters{T},
    max_upper_bound,
    step;
    kwargs...,
) where {T}
    # iterate in reverse with a negative step to find lower bound
    for r = max_upper_bound:(-step):1
        en = CircularOrbits.energy(m, r; kwargs...)
        if abs(en) > 1
            return r, max_upper_bound
        end
    end
    # for type stability
    return T(0), T(0)
end

"""
    $(TYPEDSIGNATURES)

Photon orbit radius, defined as the radius for which
```math
    \\frac{E}{\\mu} \\rightarrow \\infty .
```
"""
photon_orbit(m::AbstractMetricParameters) = error("Not implemented for $(typeof(m)).")

"""
    $(TYPEDSIGNATURES)

Marginally bound orbit
```math
    \\frac{E}{\\mu} = 1.
```
"""
marginally_bound_orbit(m::AbstractMetricParameters) =
    error("Not implemented for $(typeof(m)).")

"""
    event_horizon(m::AbstractMetricParameters; select = last, resolution = 100, θε = 1e-7, rmax = 5.0)

Event horizon radius, often equivalent to [`GradusBase.inner_radius`](@ref), however remains
distinct, such that the latter may still be an arbitrary chart cutoff.

Utility function for helping plot an event horizon shape. Returns a tuple containing the `r`
and `θ` vectors that solve

```math
    g_{t\\phi}^2 - g_{tt} g_{\\phi \\phi} = 0.
```

A `NaN` value in the `r` vector indicates no solution for that particular ``\\theta``, i.e.
that the metric describes a naked singularity.

Often the equation will have multiple roots, in which case the keyword argument `select` may be
assigned to select the desired root.
"""
event_horizon(m::AbstractMetricParameters; kwargs...) =
    error("Not implemented for $(typeof(m)).")

function _event_horizon_condition(m, r, θ)
    g = metric_components(m, (r, θ))
    g[5]^2 - g[1] * g[4]
    # inv(g[2])
end

function _solve_radius_condition(
    m::AbstractStaticAxisSymmetricParameters{T},
    condition_function;
    select = maximum,
    resolution::Int = 100,
    θε::T = T(1e-7),
    rmax = 5.0,
    init = 0.0,
) where {T}
    θs = range(θε, 2π - θε, resolution)
    rs = map(θs) do θ
        f(r) = condition_function(m, r, θ)
        r = Roots.find_zeros(f, init, rmax)
        if isempty(r) || all(isnan, r)
            NaN
        else
            select(r)
        end
    end
    rs, θs
end

function event_horizon(m::AbstractStaticAxisSymmetricParameters; kwargs...)
    _solve_radius_condition(m, _event_horizon_condition; kwargs...)
end

function is_naked_singularity(
    m::AbstractStaticAxisSymmetricParameters{T};
    resolution::Int = 100,
    θε::T = T(1e-7),
    rmax = 5.0,
) where {T}
    any(range(θε, 2π - θε, resolution)) do θ
        f(r) = _event_horizon_condition(m, r, θ)
        r = Roots.find_zeros(f, 0.0, rmax)
        isempty(r)
    end
end

function _ergosphere_condition(m, r, θ)
    g = metric_components(m, (r, θ))
    g[1]
end

function ergosphere(m::AbstractStaticAxisSymmetricParameters; kwargs...)
    _solve_radius_condition(m, _ergosphere_condition; init = 1.0, kwargs...)
end

export event_horizon, is_naked_singularity
