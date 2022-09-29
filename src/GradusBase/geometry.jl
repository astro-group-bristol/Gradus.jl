function vector_to_local_sky(m::AbstractMetricParams{T}, u, θ, ϕ) where {T}
    error("Not implemented for $(typeof(m))")
end

mdot(g, v1, v2) = @tullio r := g[i, j] * v1[i] * v2[j]
mnorm(g, v) = mdot(g, v, v)

"""
    mproject(g, v, u)

Project vector `v` onto `u` with metric `g`. Optional first argument may be
[`AbstractMetricParams`](@ref) for more optimized methods, which fallback to an einsum.
"""
mproject(g, v, u) = mdot(g, v, u) / mnorm(g, u)

function projectbasis(g, basis, v)
    s = zero(SVector{4,Float64})
    for e in basis
        s += mproject(g, v, e) .* e
    end
    s
end

function grammschmidt(v, basis, g; tol = 4eps(Float64))
    p = projectbasis(g, basis, v)

    while sum(p) > tol
        v = v .- p
        p = projectbasis(g, basis, v)
    end

    v = v .- p
    vnorm = √mnorm(g, v)
    v / vnorm
end

function tetradframe(m::AbstractMetricParams{T}, u, v) where {T}
    g = metric(m, u)

    # start procedure with ϕ, which has zero for r and θ
    vϕ = grammschmidt(@SVector[1.0, 0.0, 0.0, 1.0], (v,), g)
    # then do r, which has zero for θ
    vr = grammschmidt(@SVector[1.0, 1.0, 0.0, 1.0], (v, vϕ), g)
    # then finally θ
    vθ = grammschmidt(@SVector[1.0, 1.0, 1.0, 1.0], (v, vϕ, vr), g)

    (v, vr, vθ, vϕ)
end
