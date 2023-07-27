module __JohannsenAD
using ..StaticArrays

@fastmath f(M, r, ϵ3) = ϵ3 * M^3 / r
@fastmath Σ(r, a, θ, M, ϵ3) = r^2 + a^2 * cos(θ)^2 + f(M, r, ϵ3)
@fastmath Δ(r, M, a) = r^2 - 2M * r + a^2

@fastmath A₁(M, r, α13) = 1 + α13 * (M / r)^3
@fastmath A₂(M, r, α22) = 1 + α22 * (M / r)^2
@fastmath A₅(M, r, α52) = 1 + α52 * (M / r)^2

@fastmath function metric_components(m, rθ)
    (r, θ) = rθ

    A1 = A₁(m.M, r, m.α13)
    A2 = A₂(m.M, r, m.α22)
    A5 = A₅(m.M, r, m.α52)
    Σ₀ = Σ(r, m.a, θ, m.M, m.ϵ3)
    Δ₀ = Δ(r, m.M, m.a)

    # cache common components
    r2a2 = r^2 + m.a^2
    sin_theta2 = sin(θ)^2

    denom = ((r2a2) * A1 - m.a^2 * A2 * sin_theta2)^2

    tt = -Σ₀ * (Δ₀ - m.a^2 * A2^2 * sin_theta2)
    rr = Σ₀ / (Δ₀ * A5)
    θθ = Σ₀
    ϕϕ = Σ₀ * sin_theta2 * ((r2a2)^2 * A1^2 - m.a^2 * Δ₀ * sin_theta2)

    tϕ = -m.a * Σ₀ * sin_theta2 * ((r2a2) * A1 * A2 - Δ₀)
    @SVector [tt / denom, rr, θθ, ϕϕ / denom, tϕ / denom]
end

end # module

"""
    struct JohannsenMetric{T} <: AbstractStaticAxisSymmetric{T}

The Johannsen (20xx) metric.
- `M = 1.0`: Black hole mass.
- `a = 0.0`: Black hole spin.
- `α13 = 0.0`: ``\\alpha_{13}`` deviation parameter.
- `α22 = 0.0`: ``\\alpha_{22}`` deviation parameter.
- `α52 = 0.0`: ``\\alpha_{52}`` deviation parameter.
- `ϵ3 = 0.0`: ``\\epsilon_{3}`` deviation parameter.
"""
@with_kw struct JohannsenMetric{T} <: AbstractStaticAxisSymmetric{T}
    @deftype T
    "Black hole mass."
    M = 1.0
    "Black hole spin."
    a = 0.0
    "``\\alpha_{13}`` deviation parameter."
    α13 = 0.0
    "``\\alpha_{22}`` deviation parameter."
    α22 = 0.0
    "``\\alpha_{52}`` deviation parameter."
    α52 = 0.0
    "``\\epsilon_{3}`` deviation parameter."
    ϵ3 = 0.0
end

metric_components(m::JohannsenMetric{T}, rθ) where {T} =
    __JohannsenAD.metric_components(m, rθ)
inner_radius(m::JohannsenMetric{T}) where {T} = m.M + √(m.M^2 - m.a^2)

export JohannsenMetric
