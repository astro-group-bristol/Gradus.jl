module __BumblebeeAD
using ..StaticArrays
using ..MuladdMacro

@muladd @fastmath begin
    Δ(r, M, l) = (r^2 - (2M * r)) / (l + 1)

    function metric_components(M, a, l, rθ)
        (r, θ) = rθ
        sinθ2 = sin(θ)^2
        Δ₀ = Δ(r, M, l)

        tt = -(1 - (2M / r))
        rr = r^2 / Δ₀
        θθ = r^2
        ϕϕ = r^2 * sinθ2

        tϕ = -2M * a * sinθ2 / r
        @SVector [tt, rr, θθ, ϕϕ, tϕ]
    end
end

end # module

struct BumblebeeMetric{T} <: AbstractStaticAxisSymmetric{T}
    "Black hole mass."
    M::T
    "Black hole spin."
    a::T
    "LSB parameter."
    l::T
    function BumblebeeMetric(M::T, a, l) where {T}
        if l <= -1.0
            error("l must be >-1")
        end
        if abs(a) > 0.3
            error(
                "This metric is for the slow rotation approximation only, and requires |a| < 0.3.",
            )
        end
        new{T}(T(M), T(a), T(l))
    end
end
BumblebeeMetric(; M = 1.0, a = 0.0, l = 0.0) = BumblebeeMetric(M, a, l)

# implementation
metric_components(m::BumblebeeMetric{T}, rθ) where {T} =
    __BumblebeeAD.metric_components(m.M, m.a, m.l, rθ)
inner_radius(m::BumblebeeMetric{T}) where {T} = m.M + √(m.M^2 - m.a^2)


export BumblebeeMetric
