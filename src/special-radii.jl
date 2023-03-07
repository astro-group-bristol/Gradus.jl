"""
    $(TYPEDSIGNATURES)

Innermost stable circular orbit (ISCO), defined by
```math
    \\frac{\\text{d}}{\\text{d}r} \\left( \\frac{E}{\\mu} \\right) = 0.
```
Uses analytic solutions if known for that metric, else uses a root finder to calculate
the radius at which the above condition is met.
"""
isco(m::AbstractMetricParams) = error("Not implemented for $(typeof(m)).")

function isco(
    m::AbstractAutoDiffStaticAxisSymmetricParams,
    lower_bound,
    upper_bound;
    kwargs...,
)
    dE(r) = ForwardDiff.derivative(x -> CircularOrbits.energy(m, x; kwargs...), r)
    d2E(r) = ForwardDiff.derivative(dE, r)

    Roots.find_zero((dE, d2E), (lower_bound, upper_bound))
end

isco(
    m::AbstractAutoDiffStaticAxisSymmetricParams;
    upper_bound = 100.0,
    step = -0.005,
    kwargs...,
) = isco(
    m,
    find_lower_isco_bound(m; upper_bound = upper_bound, step = step),
    100.0;
    kwargs...,
)

"""
    $(TYPEDSIGNATURES)

Brute force method for finding the first radius at which ``\\frac{E}{\\mu} > 1``. The method
calculates ``E`` with [`CircularOrbits.energy`](@ref), where `r` steps from `upper_bound` to
``r = r_\\text{g}`` in steps of `step`.

Returns `T(0.0)` if no such radius found.
"""
function find_lower_isco_bound(
    m::AbstractAutoDiffStaticAxisSymmetricParams{T};
    upper_bound = 100.0,
    step = -0.005,
) where {T}
    # iterate in reverse with a negative step
    for r = upper_bound:step:1.0
        if CircularOrbits.energy(m, r) > 1.0
            return r
        end
    end
    # for type stability
    return T(0.0)
end

"""
    $(TYPEDSIGNATURES)

Photon orbit radius, defined as the radius for which
```math
    \\frac{E}{\\mu} \\rightarrow \\infty .
```
"""
photon_orbit(m::AbstractMetricParams) = error("Not implemented for $(typeof(m)).")

"""
    $(TYPEDSIGNATURES)

Marginally bound orbit
```math
    \\frac{E}{\\mu} = 1.
```
"""
marginally_bound_orbit(m::AbstractMetricParams) = error("Not implemented for $(typeof(m)).")

"""
    event_horizon(m::AbstractMetricParams; select = last, resolution = 100, θε = 1e-7, rmax = 5.0)

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
event_horizon(m::AbstractMetricParams; kwargs...) =
    error("Not implemented for $(typeof(m)).")

function _event_horizon_condition(m, r, θ)
    g = metric_components(m, (r, θ))
    g[5]^2 - g[1] * g[4]
    # inv(g[2])
end

function event_horizon(
    m::AbstractAutoDiffStaticAxisSymmetricParams{T};
    select = maximum,
    resolution::Int = 100,
    θε::T = T(1e-7),
    rmax = 5.0,
) where {T}
    θs = range(θε, 2π - θε, resolution)
    rs = map(θs) do θ
        f(r) = _event_horizon_condition(m, r, θ)
        r = Roots.find_zeros(f, 0.0, rmax)
        if isempty(r) || all(isnan, r)
            NaN
        else
            select(r)
        end
    end
    rs, θs
end

function is_naked_singularity(
    m::AbstractAutoDiffStaticAxisSymmetricParams{T};
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

export event_horizon, is_naked_singularity
