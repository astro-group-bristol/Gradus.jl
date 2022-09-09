module __JohannsenPsaltisAD
using ..StaticArrays

@inline h(r, M, ϵ3, Σ) = ϵ3 * M^3 * r / Σ^2
@inline Σ(r, a, θ) = r^2 + a^2 * cos(θ)^2
@inline Δ(r, M, a) = r^2 - 2M * r + a^2

@fastmath function metric_components(M, a, ϵ3, rθ)
    (r, θ) = rθ

    Σ₀ = Σ(r, a, θ)
    h₀ = h(r, M, ϵ3, Σ₀)
    sinθ2 = sin(θ)^2

    tt = -(1 + h₀) * (1 - 2M * r / Σ₀)
    rr = Σ₀ * (1 + h₀) / (Δ(r, M, a) + a^2 * sinθ2 * h₀)
    θθ = Σ₀

    term1 = sinθ2 * (r^2 + a^2 + 2a^2 * M * r * sinθ2 / Σ₀)
    term2 = h₀ * a^2 * (Σ₀ + 2M * r) * sinθ2^2 / Σ₀
    ϕϕ = term1 + term2

    tϕ = -2a * M * r * sinθ2 * (1 + h₀) / Σ₀

    @SVector [tt, rr, θθ, ϕϕ, tϕ]
end

end # module

"""
    struct JohannsenPsaltisAD{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}

Johannsen and Psaltis 2011
$(FIELDS)
"""
@with_kw struct JohannsenPsaltisAD{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}
    @deftype T
    "Black hole mass."
    M = 1.0
    "Black hole spin."
    a = 0.0
    "``\\epsilon_{3}`` deviation parameter."
    ϵ3 = 0.0
end

GeodesicTracer.metric_components(m::JohannsenPsaltisAD{T}, rθ) where {T} =
    __JohannsenPsaltisAD.metric_components(m.M, m.a, m.ϵ3, rθ)
GradusBase.inner_radius(m::JohannsenPsaltisAD{T}) where {T} = m.M + √(m.M^2 - m.a^2)