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

@fastmath function distance_to_disc(d::ThinDisc{T}, x4; gtol) where {T}
    ρ = _equatorial_project(x4)
    if ρ < d.inner_radius || ρ > d.outer_radius
        return 1.0
    end
    p = @inbounds let r = x4[2], θ = x4[3], ϕ = x4[4]
        sinθ = sin(θ)
        @SVector [r * sinθ * cos(ϕ), r * sinθ * sin(ϕ), r * cos(θ)]
    end
    n = SVector{3,T}(0, 0, 1)
    # project u into normal vector n
    k = p ⋅ n
    abs(k) - (gtol * x4[2])
end



export ThinDisc
