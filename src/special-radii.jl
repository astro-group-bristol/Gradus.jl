"""
    $(TYPEDSIGNATURES)

Innermost stable circular orbit.
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
    return 0.0
end

isco(m::AbstractAutoDiffStaticAxisSymmetricParams{T}) where {T} =
    isco(m, find_lower_isco_bound(m), 100.0)

"""
    $(TYPEDSIGNATURES)

Photon orbit.
"""
r_ph(m::AbstractMetricParams{T}) where {T} = error("Not implemented for $(typeof(m)).")

"""
    $(TYPEDSIGNATURES)

Marginally bound orbit.
"""
r_mb(m::AbstractMetricParams{T}) where {T} = error("Not implemented for $(typeof(m)).")

"""
    $(TYPEDSIGNATURES)

Event horizon.
"""
r_s(m::AbstractMetricParams{T}) where {T} = error("Not implemented for $(typeof(m)).")
