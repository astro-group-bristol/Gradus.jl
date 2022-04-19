
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

export GeometricThinDisc
