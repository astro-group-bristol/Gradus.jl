module __KerrNewmanAD
using ..StaticArrays
using ..MuladdMacro

@muladd @fastmath begin
    rQ(Q) = Q^2 / (4 * π)
    Σ(r, a, θ) = r^2 + (a * cos(θ))^2
    Δ(r, R, a, Q) = r^2 - R * r + a^2 + rQ(Q)^2

    # the way this function must be defined is a little complex
    # but helps with type-stability
    function metric_components(M, a, Q, rθ)
        (r, θ) = rθ
        R = 2M
        Σ₀ = Σ(r, a, θ)
        sinθ2 = sin(θ)^2

        tt = -(1 - (R * r) / Σ₀)
        rr = Σ₀ / Δ(r, R, a, Q)
        θθ = Σ₀
        ϕϕ = sinθ2 * (r^2 + a^2 + (sinθ2 * R * r * a^2) / Σ₀)

        tϕ = (-R * r * a * sinθ2) / Σ₀
        @SVector [tt, rr, θθ, ϕϕ, tϕ]
    end

    function electromagnetic_potential(m, rθ)
        (r, θ) = rθ
        Σ₀ = Σ(r, m.a, θ)
        rq = rQ(m.Q)
        SVector{4,eltype(rθ)}(r * rq / Σ₀, 0, 0, -m.a * r * rq * sin(θ)^2 / Σ₀)
    end
end

end # module

"""
    struct KerrNewmanMetric{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}
"""
struct KerrNewmanMetric{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}
    "Black hole mass."
    M::T
    "Black hole spin."
    a::T
    "Black hole charge."
    Q::T

    function KerrNewmanMetric(M::T, a, Q) where {T}
        if a^2 + Q^2 > M^2
            error("Value error: `a^2 + Q^2` must be `<= M^2`")
        end
        new{T}(T(M), T(a), T(Q))
    end
end

# implementation
metric_components(m::KerrNewmanMetric{T}, rθ) where {T} =
    __KerrNewmanAD.metric_components(m.M, m.a, m.Q, rθ)
inner_radius(m::KerrNewmanMetric{T}) where {T} =
    m.M + √(m.M^2 - m.a^2 - __KerrNewmanAD.rQ(m.Q)^2)

electromagnetic_potential(m::KerrNewmanMetric, x) =
    __KerrNewmanAD.electromagnetic_potential(m, x)

function geodesic_ode_problem(
    m::KerrNewmanMetric{T},
    pos::StaticVector{S,T},
    vel::StaticVector{S,T},
    time_domain,
    ;
    q = 0.0,
    kwargs...,
) where {S,T}
    function f(u::SVector{8,T}, p, λ) where {T}
        @inbounds let x = SVector{4,T}(@view(u[1:4])), v = SVector{4,T}(@view(u[5:8]))
            dv = SVector{4,T}(geodesic_eq(m, x, v))
            # add maxwell part
            F = maxwell_tensor(m, x)
            dvf = q * (F * v)
            vcat(v, dv + dvf)
        end
    end

    u_init = vcat(pos, vel)
    ODEProblem{false}(
        f,
        u_init,
        time_domain,
        IntegrationParameters(StatusCodes.NoStatus);
        kwargs...,
    )
end

export KerrNewmanMetric
