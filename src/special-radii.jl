"""
    $(TYPEDSIGNATURES)

Innermost stable circular orbit.
"""
isco(m::AbstractMetricParams{T}) where {T} = error("Not implemented for $(typeof(m)).")

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
