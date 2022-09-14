"""
    intersects_geometry(m::AbstractAccretionGeometry{T}, line_element)

Utility function. Returns a boolean dependent on whether `line_element` intersects with the geometry (`true`) or not (`false`).
Uses [`in_nearby_region`](@ref) to optimze and calls [`has_intersect`](@ref) to determine intersection.
"""
function intersects_geometry(
    m::AbstractAccretionGeometry{T},
    line_element,
    integrator,
) where {T}
    if in_nearby_region(m, line_element)
        return has_intersect(m, line_element)
    end
    false
end

"""
    in_nearby_region(g::AbstractAccretionGeometry{T}, line_element)::Bool

Returns a boolean indicating whether the `line_element` is "close" to the geometry in
question. Used to optimize when to call the intersection algorithm.

Contextually, "close" is a little arbitrary, and this function may always return `True`, however
will suffer in performance if this is the case.

`line_element` is a tuple of two four-position vectors, indicating the last integration position,
and current integration position, i.e. `(u_prev, u_curr)`, in integrator coordinates.

# Notes

This function actually depends on the step size of the integrator, but this is currently not
considered in the implementation.
"""
in_nearby_region(d::AbstractAccretionGeometry{T}, line_element) where {T} =
    error("Not implemented for $(typeof(d))")

"""
    has_intersect(g::AbstractAccretionGeometry{T}, line_element)

Returns a boolean indicating whether `line_element` intersects the geometry `g`. The intersection
algorithm used depends on the geometry considered. For meshes, this uses the [`jsf_algorithm`](@ref).

`line_element` is a tuple of two four-position vectors, indicating the last integration position,
and current integration position, i.e. `(u_prev, u_curr)`, in integrator coordinates.
"""
has_intersect(d::AbstractAccretionGeometry{T}, line_element) where {T} =
    error("Not implemented for $(typeof(d))")

"""
    jsf_algorithm(V₁::T, V₂::T, V₃::T, Q₁::V, Q₂::V; ϵ = 1e-8)

Implemented from Jiménez, Segura, Feito. Computation Geometry 43 (2010) 474-492.

See [this blog post](https://fjebaker.github.io/blog/pages/2022-01-ray-tracing-a-cow/#jim%C3%A9nez_segura_and_feito_2010)
for a discussion.
"""
function jsf_algorithm(V₁::T, V₂::T, V₃::T, Q₁::V, Q₂::V; ϵ = 1e-8) where {T,V}
    A = Q₁ .- V₃
    B = V₁ .- V₃
    C = V₂ .- V₃
    W₁ = B × C
    w = A ⋅ W₁
    if w > ϵ
        D = Q₂ .- V₃
        s = D ⋅ W₁
        s > ϵ && return false, float(0)
        W₂ = A × D
        t = W₂ ⋅ C
        t < -ϵ && return false, float(0)
        u = -W₂ ⋅ B
        u < -ϵ && return false, float(0)
        w < s + t + u && return false, float(0)
    elseif w < -ϵ
        return false, float(0)
    else # w == 0
        D = Q₂ .- V₃
        s = D ⋅ W₁
        if s > ϵ
            return false, float(0)
        elseif s < -ϵ
            W₂ = D × A
            t = W₂ ⋅ C
            t > ϵ && return false, float(0)
            u = -W₂ ⋅ B
            u > ϵ && return false, float(0)
            -s > t + u && return false, float(0)
        else
            return false, float(0)
        end
    end
    t = w / (s - w)
    return true, float(t)
end
