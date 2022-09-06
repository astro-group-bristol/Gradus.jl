function vector_to_local_sky(m::AbstractMetricParams{T}, u, θ, ϕ) where {T}
    error("Not implemented for $(typeof(m))")
end

mdot(g, v1, v2) = @tullio r := g[i, j] * v1[i] * v2[j]
mnorm(g, v) = mdot(g, v, v)

# fallback methods
mdot(::AbstractMetricParams{T}, g, v1, v2) where {T} = mdot(g, v1, v2)
mnorm(m::AbstractMetricParams{T}, g, v) where {T} = mdot(m, g, v, v)

"""
    mproject(g, v, u)
    mproject(m::AbstractMetricParams{T}, g, v, u)

Project vector `v` onto `u` with metric `g`. Optional first argument may be
[`AbstractMetricParams`](@ref) for more optimized methods, which fallback to an einsum.
"""
mproject(m::AbstractMetricParams{T}, g, v, u) where {T} = mdot(m, g, v, u) / mnorm(m, g, u)
mproject(g, v, u) = mdot(g, v, u) / mnorm(g, u)

function projectbasis(m::AbstractMetricParams{T}, g, basis, v::AbstractArray{T}) where {T}
    s = zero(SVector{4,Float64})
    for e in basis
        s += mproject(m, g, v, e) .* e
    end
    s
end

function grammschmidt(m::AbstractMetricParams{T}, g, basis; tol = eps(T)) where {T}
    v = ones(SVector{4,Float64})
    p = projectbasis(m, g, basis, v)

    while sum(p) > 4 * tol
        v = v .- p
        p = projectbasis(m, g, basis, v)
    end

    v = v .- p
    v ./ √mnorm(m, g, v)
end

function tetradframe(m::AbstractMetricParams{T}, u, v) where {T}
    g = metric(m, u)
    M = zero(MMatrix{4,4,T})
    known_frame_vectors = minimalframe(m, g, u, v)
    for (i, e) in enumerate(known_frame_vectors)
        M[:, i] .= e
    end
    for i = length(known_frame_vectors)+1:4
        M[:, i] .= grammschmidt(m, g, eachcol(M[:, 1:i-1]))
    end
    M
end

function minimalframe(m::AbstractMetricParams{T}, g, u, v) where {T}
    (v,)
end
