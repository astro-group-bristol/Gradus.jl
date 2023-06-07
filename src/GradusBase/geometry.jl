function vector_to_local_sky(m::AbstractMetric, u, θ, ϕ)
    error("Not implemented for $(typeof(m))")
end

dotproduct(g::AbstractMatrix, v1, v2) = @tullio r := g[i, j] * v1[i] * v2[j]
propernorm(g::AbstractMatrix, v) = dotproduct(g, v, v)
propernorm(m::AbstractMetric, u, v) = propernorm(metric(m, u), v)

"""
    mproject(g, v, u)

Project vector `v` onto `u` with metric `g`. Optional first argument may be
[`AbstractMetric`](@ref) for more optimized methods, which fallback to an einsum.
"""
mproject(g, v, u) = dotproduct(g, v, u) / propernorm(g, u)

function projectbasis(g, basis, v)
    s = zero(SVector{4,eltype(g)})
    for e in basis
        s += mproject(g, v, e) .* e
    end
    s
end

function gramschmidt(v, basis, g; tol = 4eps(eltype(g)))
    p = projectbasis(g, basis, v)

    while sum(p) > tol
        v = v .- p
        p = projectbasis(g, basis, v)
    end

    v = v .- p
    vnorm = √abs(propernorm(g, v))
    v / vnorm
end

# TODO: this presupposes static and axis symmetric
@inline function tetradframe(g::AbstractMatrix{T}, v) where {T}
    vt = v ./ √abs(propernorm(g, v))
    # start procedure with ϕ, which has zero for r and θ
    vϕ = gramschmidt(SVector{4,T}(1, 0, 0, 1), (vt,), g)
    # then do r, which has zero for θ
    vr = gramschmidt(SVector{4,T}(0, 1, 0, 0), (vt, vϕ), g)
    # then finally θ
    vθ = gramschmidt(SVector{4,T}(0, 0, 1, 0), (vt, vϕ, vr), g)
    (vt, vr, vθ, vϕ)
end

tetradframe(m::AbstractMetric, u, v) = tetradframe(metric(m, u), v)

# TODO: this presupposes static and axis symmetric
# tetrad with latin indices down: frame
function lnrframe(g::AbstractMatrix{T}) where {T}
    ω = -g[1, 4] / g[4, 4]
    v = SVector{4,T}(1, 0, 0, ω)
    tetradframe(g, v)
end
lnrframe(m::AbstractMetric, u) = lnrframe(metric(m, u))

# tetrad with latin indices up: basis
function lnrbasis(g::AbstractMatrix{T}) where {T}
    ω = -g[1, 4] / g[4, 4]
    v = SVector{4,T}(-ω, 0, 0, 1)
    (vϕ, vr, vθ, vt) = tetradframe(inv(g), v)
    # rearrange
    (vt, vr, vθ, vϕ)
end
lnrbasis(m::AbstractMetric, u) = lnrbasis(metric(m, u))

lowerindices(g::AbstractMatrix, v) = g * v
lowerindices(m::AbstractMetric, u, v) = lowerindices(metric(m, u), v)

raiseindices(ginv::AbstractMatrix, v) = ginv * v
raiseindices(m::AbstractMetric, u, v) = raiseindices(inv(metric(m, u)), v)
