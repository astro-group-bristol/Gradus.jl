"""
Einstein-Maxwell-Dilaton-Axion (EMDA) metric.
García et al. (1995).
"""
module __DilatonAxionAD
using ..StaticArrays

Σ(r, a, θ) = r^2 + a^2 * cos(θ)^2
Δ(r, R, a) = r^2 + a^2 - R * r
# the a here is actually the dimensionless spin parameter a / R ?
Σhat(Σ, β, b, r, βb, a, θ, R) = Σ - (β^2 + 2b * r) + R^2 * βb * (βb - 2 * a * cos(θ))
Δhat(Δ, β, b, r, βb, R) = Δ - (β^2 + 2b * r) - R * (R + 2b) * βb^2

δ(r, b, a) = r^2 - 2b * r + a^2

W(θ, βab, βa) = 1 + (βab * (2 * cos(θ) - βab) + βa^2) * csc(θ)^2
A(a, θ, δ, W, Δhat) = δ^2 - Δhat * (W * a * sin(θ))^2

@fastmath function metric_components(M, a, β, b, rθ)
    (r, θ) = rθ

    R = M

    # check these for the normalisation term R
    βb = β / b
    βa = β / a
    βab = β / (a * b)

    Σ₀ = Σ(r, a, θ)
    Δ₀ = Δ(r, R, a)
    Δhat₀ = Δhat(Δ₀, β, b, r, βb, R)
    Σhat₀ = Σhat(Σ₀, β, b, r, βb, a, θ, R)
    δ₀ = δ(r, b, a)
    W₀ = W(θ, βab, βa)
    A₀ = A(a, θ, δ₀, W₀, Δhat₀)

    tt = -(Δhat₀ - a^2 * sin(θ)^2) / Σhat₀
    rr = Σhat₀ / Δhat₀
    θθ = Σhat₀
    ϕϕ = A₀ * sin(θ)^2 / Σhat₀

    tϕ = -a * (δ₀ - Δhat₀ * W₀) * sin(θ)^2 / Σhat₀
    @SVector [tt, rr, θθ, ϕϕ, tϕ]
end

end # module

"""
    struct DilatonAxion{T} <: AbstractStaticAxisSymmetric{T}

Einstein-Maxwell-Dilaton-Axion metric.
- `M = 1.0`: Singularity mass.
- `a = 0.0`: Singularity spin.
- `β = 0.0`: Dilaton coupling strength.
- `b = 1.0`: Axion coupling strength.
"""
@with_kw struct DilatonAxion{T} <: AbstractStaticAxisSymmetric{T}
    @deftype T
    "Singularity mass."
    M = 1.0
    "Singularity spin."
    a = 0.5
    "Dilaton coupling strength."
    β = 0.0
    "Axion coupling strength."
    b = 1.0
end

metric_components(m::DilatonAxion{T}, rθ) where {T} =
    __DilatonAxionAD.metric_components(m.M, m.a, m.β, m.b, rθ)
inner_radius(m::DilatonAxion{T}) where {T} =
    m.M + m.b + √((m.M + m.b)^2 - m.a^2 + m.β^2 - (m.M - 2m.b) * m.M * (m.β / m.b)^2)

export DilatonAxion
