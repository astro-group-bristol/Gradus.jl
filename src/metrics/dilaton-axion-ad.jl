"""
Einstein-Maxwell-Dilaton-Axion (EMDA) metric.
García et al. (1995).
"""
module __DilatonAxionAD
using ..StaticArrays

@fastmath Σ(r, a, θ) = r^2 + a^2 * cos(θ)^2
@fastmath Δ(r, a, rg) = r^2 - 2rg * r + a^2
@fastmath A(r, a, Δ, θ) = (r^2 + a^2)^2 - a^2 * Δ * sin(θ)^2

@fastmath W(βab, θ, βa) = 1 + (βab * (2 * cos(θ) - βab) + βa^2) * csc(θ)^2
@fastmath Σ̂(Σ, β, b, r, βb, a, θ, rg) =
    Σ - (β^2 + 2b * r) + rg^2 * βb * (βb - 2 * (a / rg) * cos(θ))
@fastmath Δ̂(Δ, β, b, r, βb, rg) = Δ - (β^2 + 2b * r) - rg * (rg + 2b) * βb^2
@fastmath Â(δ, a, Δ̂, W, θ) = δ^2 - a^2 * Δ̂ * W^2 * sin(θ)^2
@fastmath δ(r, b, a) = r^2 - 2b * r + a^2

@fastmath function metric_components(M, a, β, b, rθ)
    (r, θ) = rθ

    rg = M

    Δ₀ = Δ(r, a, rg)
    Σ₀ = Σ(r, a, θ)

    βa = β / a
    βb = β / b
    βab = rg * β / (a * b)

    Δ̂₀ = Δ̂(Δ₀, β, b, r, βb, rg)
    Σ̂₀ = Σ̂(Σ₀, β, b, r, βb, a, θ, rg)
    δ₀ = δ(r, b, a)
    W₀ = W(βab, θ, βa)

    Â₀ = Â(δ₀, a, Δ̂₀, W₀, θ)

    tt = -(Δ̂₀ - a^2 * sin(θ)^2) / Σ̂₀
    rr = Σ̂₀ / Δ̂₀
    θθ = Σ̂₀
    ϕϕ = Â₀ * sin(θ)^2 / Σ̂₀

    tϕ = -a * (δ₀ - Δ̂₀ * W₀) * sin(θ)^2 / Σ̂₀
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
    a = 0.0
    "Dilaton coupling strength."
    β = 0.0
    "Axion coupling strength."
    b = 1.0
end

metric_components(m::DilatonAxion{T}, rθ) where {T} =
    __DilatonAxionAD.metric_components(m.M, m.a, m.β, m.b, rθ)
inner_radius(m::DilatonAxion{T}) where {T} = m.M + √(m.M^2 - m.a^2)

export DilatonAxion
