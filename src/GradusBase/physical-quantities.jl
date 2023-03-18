
# both energy and angular momentum
# assume time only coupled to radial coordinate
# need to think of a nice way to keep this efficient
# whilst allowing metrics with other couplings

"""
    E(m::AbstractMatrix{T}, v)
    E(m::AbstractMetricParameters{T}, u, v)

Compute the energy for a numerically evaluated metric, and some velocity four vector `v`,
```math
E = - p_t = - g_{t\\nu} p^\\nu.
```

For null geodesics, the velocity is the momentum ``v^\\nu = p^\\nu``. For massive geodesics,
the mass ``\\mu`` needs to be known to compute ``\\mu v^\\nu = p^\\nu``.
"""
function E(metric::AbstractMatrix{T}, v) where {T}
    T(@inbounds -(metric[1, 1] * v[1] + metric[1, 4] * v[4]))
end
E(m::AbstractMetricParameters, u, v) = E(metric(m, u), v)
E(m::AbstractMetricParameters, gp::AbstractGeodesicPoint) = E(m, gp.u2, gp.v2)


"""
    Lz(m::AbstractMatrix{T}, v)
    Lz(m::AbstractMetricParameters{T}, u, v)

Compute the angular momentum for a numerically evaluated metric, and some velocity four vector `v`.
```math
L_z = p_\\phi = - g_{\\phi\\nu} p^\\nu.
```
"""
function Lz(metric::AbstractMatrix{T}, v) where {T}
    T(@inbounds metric[4, 4] * v[4] + metric[1, 4] * v[1])
end
Lz(m::AbstractMetricParameters, u, v) = Lz(metric(m, u), v)
Lz(m::AbstractMetricParameters, gp::AbstractGeodesicPoint) = Lz(m, gp.u2, gp.v2)
