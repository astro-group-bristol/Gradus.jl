module __KerrRefractiveAD
using ..StaticArrays

@fastmath Σ(r, a, θ) = r^2 + a^2 * cos(θ)^2
@fastmath Δ(r, R, a) = r^2 - R * r + a^2

@fastmath function metric_components(M, a, n, corona_radius, rθ)
    (r, θ) = rθ
    R = 2M
    Σ₀ = Σ(r, a, θ)
    sinθ2 = sin(θ)^2

    Δ₀ = Δ(r, R, a)

    tt = -(1 - (R * r) / Σ₀)
    rr = Σ₀ / Δ₀
    θθ = Σ₀
    ϕϕ = sinθ2 * (r^2 + a^2 + (sinθ2 * R * r * a^2) / Σ₀)

    tϕ = (-R * r * a * sinθ2) / Σ₀

    # need to interpolate the refractve index boundary
    # else we don't have a gradient to compute, and we miss it
    # hemorraging energy in the process
    δr = 2.0
    if r ≤ corona_radius
        # ...
    elseif corona_radius ≤ r ≤ corona_radius + δr
        t = (r - corona_radius) / δr
        # use an arbitrarily steep smooth interpolation
        # this one isn't perfect, but does a good job
        k = atan(1e5t) * 2 / π
        n = k + n * (1 - k)
    else
        n = 1.0
    end

    @SVector [tt / n^2, rr, θθ, ϕϕ, tϕ / n]
end

end # module

"""
    struct KerrRefractive{T} <: AbstractStaticAxisSymmetric{T}

Kerr metric in Boyer-Lindquist coordintes with a path-length ansatz, equivalent to a refractive
index `n`, within the coronal radius `corona_radius`.
- `M = 1.0`: Black hole mass.
- `a = 0.0`: Black hole spin.
- `n = 1.0`: Refractive index within the corona.
- `corona_radius = 20.0`: Radius of the corona.
"""
@with_kw struct KerrRefractive{T} <: AbstractStaticAxisSymmetric{T}
    @deftype T
    "Black hole mass."
    M = 1.0
    "Black hole spin."
    a = 0.0
    "Refractive index within the corona."
    n = 1.0
    "Radius of the corona."
    corona_radius = 20.0
end

# implementation
metric_components(m::KerrRefractive{T}, rθ) where {T} =
    __KerrRefractiveAD.metric_components(m.M, m.a, m.n, m.corona_radius, rθ)
inner_radius(m::KerrRefractive{T}) where {T} = m.M + √(m.M^2 - m.a^2)

# special radii
isco(m::KerrRefractive{T}) where {T} = __BoyerLindquistFO.isco(m.M, m.a)

export KerrRefractive
