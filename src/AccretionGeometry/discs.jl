
"""
    struct GeometricThinDisc{T} <: AbstractAccretionDisc{T}
    GeometricThinDisc(inner_radius::T, outer_radius::T, inclination::T)

$(FIELDS)

Simple geometrically thin accretion disc spanning from `inner_radius` to `outer_radius` in
gravitational units. Inclination of the disc is relative to spin axis, with ``90^\\circ``
being perpendicular to the spin axis.
"""
struct GeometricThinDisc{T} <: AbstractAccretionDisc{T}
    inner_radius::T
    outer_radius::T
    inclination::T
end

function in_nearby_region(d::GeometricThinDisc{T}, line_element) where {T}
    p = line_element[2]
    d.inner_radius < p[1] < d.outer_radius
end

function has_intersect(d::GeometricThinDisc{T}, line_element) where {T}
    u1, u2 = line_element
    sinα = sin(d.inclination)
    cosα = cos(d.inclination)

    s1 = u1[1] * (cosα * sin(u1[2]) * cos(u1[3]) + sinα * cos(u1[2]))
    s2 = u2[1] * (cosα * sin(u2[2]) * cos(u2[3]) + sinα * cos(u2[2]))

    s1 * s2 ≤ 0
end

function build_collision_callback(g::GeometricThinDisc{T}; gtol) where {T}
    ContinuousCallback(
        (u, λ, integrator) -> distance_to_disc(g, u; gtol = gtol),
        i -> terminate!(i, :Intersected);
        interp_points = 8,
        save_positions = (true, false),
    )
end

@fastmath function distance_to_disc(m::GeometricThinDisc{T}, u4; gtol) where {T}
    p = @inbounds let r = u4[2], θ = u4[3], ϕ = u4[4]
        if r < m.inner_radius || r > m.outer_radius
            return 1.0
        end
        @SVector [r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ)]
    end
    n = @SVector [T(0.0), cos(m.inclination), sin(m.inclination)]
    # project u into normal vector n
    k = p ⋅ n
    abs(k) - (gtol * u4[2])
end
