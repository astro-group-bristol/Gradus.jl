
function intersects_geometry(m::AbstractAccretionGeometry{T}, line_element) where {T}
    if in_nearby_region(m, line_element)
        return has_intersect(m, line_element)
    end
    false
end

in_nearby_region(d::AbstractAccretionGeometry{T}, line_element) where {T} =
    error("Not implemented for $(typeof(d))")
has_intersect(d::AbstractAccretionGeometry{T}, line_element) where {T} =
    error("Not implemented for $(typeof(d))")

# Jiménez, Segura, Feito. Computation Geometry 43 (2010) 474-492
function jsr_algorithm(V₁::T, V₂::T, V₃::T, Q₁::V, Q₂::V; ϵ = 1e-8) where {T,V}
    A = Q₁ .- V₃
    B = V₁ .- V₃
    C = V₂ .- V₃
    W₁ = B × C
    w = A ⋅ W₁
    if w > ϵ
        D = Q₂ .- V₃
        s = D ⋅ W₁
        s > ϵ && return false
        W₂ = A × D
        t = W₂ ⋅ C
        t < -ϵ && return false
        u = -W₂ ⋅ B
        u < -ϵ && return false
        w < s + t + u && return false
    elseif w < -ϵ
        return false
    else # w == 0
        D = Q₂ .- V₃
        s = D ⋅ W₁
        if s > ϵ
            return false
        elseif s < -ϵ
            W₂ = D × A
            t = W₂ ⋅ C
            t > ϵ && return false
            u = -W₂ ⋅ B
            u > ϵ && return false
            -s > t + u && return false
        else
            return false
        end
    end
    true
end
