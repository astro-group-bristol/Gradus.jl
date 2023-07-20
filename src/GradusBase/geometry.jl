# use this everywhere where we need a dot product so it's quick and easy to
# change the underlying implementation
@fastmath function _fast_dot(x, y)
    @assert size(x) == size(y)
    res = zero(eltype(x))
    @inbounds for i in eachindex(x)
        res += x[i] * y[i]
    end
    res
end

function vector_to_local_sky(m::AbstractMetric, u, θ, ϕ)
    error("Not implemented for $(typeof(m))")
end

dotproduct(g::AbstractMatrix, v1, v2) = _fast_dot(g * v1, v2)
dotproduct(m::AbstractMetric, x, v1, v2) = dotproduct(metric(m, x), v1, v2)
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

"""
    _tetrad_permute(x)

Permute the spatial components (vectors) of a given tetrad. That is, given
some `x = (e1, e2, e3, e4)`, a single call to `_tetrad_permute` will return
`(e1, e4, e2, e3)`.

Attempts to respect the input type in the return type. Checked for `SVector` and `NTuple`.
"""
_tetrad_permute(x::SVector{4,T}) where {T} = SVector{4,T}(x[1], x[4], x[2], x[3])
_tetrad_permute(x::NTuple{4}) = (x[1], x[4], x[2], x[3])

"""
    tetradframe(m::AbstractMetric, x, v)
    tetradframe(g::AbstractMatrix [= metric(m, x)], v)

Compute a tetrad frame via the Gram-Schmidt orthonormalization procedure (see [`gramschmidt`](@ref)),
where the first vector of the frame is the velocity vector `v` at some position `x`.

Used to compute the locally non-rotating frame via [`lnrframe`](@ref) and [`lnrbasis`](@ref) 
respectively.

Attempts to reorder the tetrad into a canonical form associated with the global coordinates,
that is, returns a tuple that corresponds to the ``x^1, x^2, x^3, x^4`` coordinates of the spacetime.
"""
@inline function tetradframe(g::AbstractMatrix{T}, v) where {T}
    # TODO: this presupposes static and axis symmetric
    # normalise the vector to unit length
    v1 = v ./ √abs(propernorm(g, v))

    state = v1 .!= 0
    # ensure there is an initial direction
    # else just set it to r
    if sum(state) == 1
        state = SVector{4,eltype(state)}(1, 0, 0, 1)
    end
    # store number of permutations of the space vectors needed to get 
    # tetrad in the correct order at the end
    permutations = searchsortedfirst(@views(state[2:end]), 1)
    v2 = gramschmidt(SVector{4,T}(state), (v1,), g)

    state = state .| _tetrad_permute(state)
    v3 = gramschmidt(SVector{4,T}(state), (v1, v2), g)

    state = state .| _tetrad_permute(state)
    v4 = gramschmidt(SVector{4,T}(state), (v1, v2, v3), g)

    # order and return t, r, θ, ϕ
    ret = (v1, v2, v3, v4)
    for _ = 2:permutations
        ret = _tetrad_permute(ret)
    end
    ret
end
tetradframe(m::AbstractMetric, x, v) = tetradframe(metric(m, x), v)

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

lnrbasis_matrix(m::AbstractMetric, x) = reduce(hcat, lnrbasis(m, x))
lnrframe_matrix(m::AbstractMetric, x) = reduce(hcat, lnrframe(m, x))

export dotproduct, propernorm, _fast_dot
