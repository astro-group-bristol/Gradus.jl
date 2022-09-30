function vector_to_local_sky(m::AbstractMetricParams, u, θ, ϕ)
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

function gramschmidt(v, basis, g; tol = 4eps(Float64))
    p = projectbasis(g, basis, v)

    while sum(p) > tol
        v = v .- p
        p = projectbasis(g, basis, v)
    end

    v = v .- p
    vnorm = √abs(mnorm(g, v))
    v / vnorm
end

# TODO: this presupposes static and axis symmetric
@inline function _tetradframe(g, v)
    vt = v ./ √abs(mnorm(g, v))
    # start procedure with ϕ, which has zero for r and θ
    vϕ = gramschmidt(@SVector[1.0, 0.0, 0.0, 1.0], (vt,), g)
    # then do r, which has zero for θ
    vr = gramschmidt(@SVector[0.0, 1.0, 0.0, 0.0], (vt, vϕ), g)
    # then finally θ
    vθ = gramschmidt(@SVector[0.0, 0.0, 1.0, 0.0], (vt, vϕ, vr), g)
    (vt, vr, vθ, vϕ)
end

tetradframe(m::AbstractMetricParams, u, v) = _tetradframe(metric(m, u), v)

# TODO: this presupposes static and axis symmetric
# tetrad with indices down: frame
function lnrframe(m::AbstractMetricParams, u)
    g = metric(m, u)
    ω = -g[1, 4] / g[4, 4]
    v = @SVector [1.0, 0.0, 0.0, ω]
    _tetradframe(g, v)
end

# tetrad with indices up: basis
function lnrbasis(m::AbstractMetricParams, u)
    g = metric(m, u)
    ω = -g[1, 4] / g[4, 4]
    v = @SVector [-ω, 0.0, 0.0, 1.0]
    (vϕ, vr, vθ, vt) = _tetradframe(inv(g), v)
    # rearrange
    (vt, vr, vθ, vϕ)
end

lowerindices(g, v) = g * v
lowerindices(m::AbstractMetricParams, u, v) = lowerindices(metric(m, u), v)

raiseindices(ginv, v) = ginv * v
raiseindices(m::AbstractMetricParams, u, v) = raiseindices(inv(metric(m, u)), v)
