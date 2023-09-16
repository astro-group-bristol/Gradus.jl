"""
    struct ThinDisc{T} <: AbstractAccretionDisc{T}
    ThinDisc(inner_radius::T, outer_radius::T)

Simple geometrically thin accretion disc spanning from `inner_radius` to `outer_radius` in
gravitational units. Inclination of the disc is relative to spin axis, with ``90^\\circ``
being perpendicular to the spin axis.
"""
struct ThinDisc{T} <: AbstractAccretionDisc{T}
    inner_radius::T
    outer_radius::T
end

ThinDisc(; inner_radius = 0.0, outer_radius = 500.0) = ThinDisc(inner_radius, outer_radius)

optical_property(::Type{<:ThinDisc}) = OpticallyThin()

inner_radius(disc::ThinDisc) = disc.inner_radius

function distance_to_disc(d::ThinDisc, x4; gtol)
    ρ = _equatorial_project(x4)
    if ρ < d.inner_radius || ρ > d.outer_radius
        return one(eltype(x4))
    end
    _spinaxis_project(x4, signed = false) - _gtol_error(gtol, x4)
end

export ThinDisc
