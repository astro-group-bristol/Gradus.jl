"""
    $(TYPEDSIGNATURES)

Innermost stable circular orbit (ISCO), defined by
```math
    \\frac{\\text{d}}{\\text{d}r} \\left( \\frac{E}{\\mu} \\right) = 0.
```
Uses analytic solutions if known for that metric, else uses a root finder to calculate
the radius at which the defining condition is met.
"""
isco(m::AbstractMetricParams{T}) where {T} = error("Not implemented for $(typeof(m)).")

function isco(
    m::AbstractAutoDiffStaticAxisSymmetricParams{T},
    lower_bound,
    upper_bound,
) where {T}
    dE(r) = ForwardDiff.derivative(x -> CircularOrbits.energy(m, x), r)
    d2E(r) = ForwardDiff.derivative(dE, r)

    Roots.find_zero((dE, d2E), (lower_bound, upper_bound))
end

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

isco(m::AbstractAutoDiffStaticAxisSymmetricParams{T}) where {T} =
    isco(m, find_lower_isco_bound(m), 100.0)

"""
    $(TYPEDSIGNATURES)

Photon orbit radius, defined as the radius for which
```math
    \\frac{E}{\\mu} \\rightarrow \\infty .
```
"""
r_ph(m::AbstractMetricParams{T}) where {T} = error("Not implemented for $(typeof(m)).")

"""
    $(TYPEDSIGNATURES)

Marginally bound orbit
```math
    \\frac{E}{\\mu} = 1.
```
"""
r_mb(m::AbstractMetricParams{T}) where {T} = error("Not implemented for $(typeof(m)).")

"""
    $(TYPEDSIGNATURES)

Event horizon radius, often equivalent to [`GradusBase.inner_radius`](@ref), however remains
distinct, such that the latter may still be an arbitrary chart cutoff.
"""
r_s(m::AbstractMetricParams{T}) where {T} = error("Not implemented for $(typeof(m)).")

"""
    eventhorizon(m::AbstractMetricParams; select = last, resolution = 100, θε = 1e-7, rmax = 5.0)

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
eventhorizon(m::AbstractMetricParams{T}; kwargs...) where {T} =
    error("Not implemented for $(typeof(m)).")

function _event_horizon_condition(m, r, θ)
    g = metric_components(m, (r, θ))
    g[5]^2 - g[1] * g[4]
end

function eventhorizon(
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
        if isempty(r)
            push!(r, NaN)
        end
        select(r)
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

export eventhorizon, is_naked_singularity
