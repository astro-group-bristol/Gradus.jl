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

"""
    struct WarpedThinDisc{T} <: AbstractAccretionDisc{T}
    WarpedThinDisc(inner_radius::T, outer_radius::T, f::F)

Similar to [`ThinDisc`](@ref), however the scale height as a function of the radial coordinate may be specified by an arbitrary function. The function should return the ``h = \\cos(\\theta)`` signed scaled height (c.f. [`ThickDisc`](@ref) which is unsigned).

The function should be of the form

```julia
function scale_height(ρ)
    # ...
end
```
"""
struct WarpedThinDisc{T,F} <: AbstractAccretionDisc{T}
    f::F
    inner_radius::T
    outer_radius::T
end

WarpedThinDisc(f; inner_radius = 0.0, outer_radius = 500.0) = WarpedThinDisc(f, inner_radius, outer_radius)

optical_property(::Type{<:WarpedThinDisc}) = OpticallyThin()

inner_radius(disc::WarpedThinDisc) = disc.inner_radius

function distance_to_disc(d::WarpedThinDisc, x4; gtol)
    ρ = _equatorial_project(x4)
    if ρ < d.inner_radius || ρ > d.outer_radius
        return one(eltype(x4))
    end

    h = d.f(ρ)
    γ = _spinaxis_project(x4, signed = true)

    abs(h - γ) - _gtol_error(gtol, x4)
end

export ThinDisc, WarpedThinDisc

