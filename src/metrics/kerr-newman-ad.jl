module __KerrNewmanAD
using ..StaticArrays
using ..MuladdMacro

@muladd @fastmath begin
    Σ(r, a, θ) = r^2 + (a * cos(θ))^2
    Δ(r, R, a, Q) = r^2 - R * r + a^2 + Q^2

    # the way this function must be defined is a little complex
    # but helps with type-stability
    function metric_components(M, a, Q, rθ)
        (r, θ) = rθ
        R = 2M
        Σ₀ = Σ(r, a, θ)
        sinθ2 = sin(θ)^2
        Δ₀ = Δ(r, R, a, Q)

        r2a2 = r^2 + a^2
        tt = (a^2 * sinθ2 - Δ₀) / Σ₀
        rr = Σ₀ / Δ₀
        θθ = Σ₀
        ϕϕ = (sinθ2 / Σ₀) * (r2a2^2 - a^2 * sinθ2 * Δ₀)

        tϕ = (a * sinθ2 / Σ₀) * (Δ₀ - r2a2)
        @SVector [tt, rr, θθ, ϕϕ, tϕ]
    end

    function electromagnetic_potential(m, rθ)
        (r, θ) = rθ
        Σ₀ = Σ(r, m.a, θ)
        (r * m.Q / Σ₀) * SVector{4,eltype(rθ)}(1, 0, 0, -m.a * sin(θ)^2)
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

KerrNewmanMetric(; M = 1.0, a = 0.0, Q = 0.0) = KerrNewmanMetric(M, a, Q)

# implementation
metric_components(m::KerrNewmanMetric{T}, rθ) where {T} =
    __KerrNewmanAD.metric_components(m.M, m.a, m.Q, rθ)
inner_radius(m::KerrNewmanMetric{T}) where {T} = m.M + √(m.M^2 - m.a^2 - m.Q^2)

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
            dvf = if !(q ≈ 0.0)
                F = faraday_tensor(m, x)
                q * (F * v)
            else
                SVector{4,T}(0, 0, 0, 0)
            end
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

function CircularOrbits.energy(m::KerrNewmanMetric, rθ, utuϕ; q = 0.0)
    V = __KerrNewmanAD.electromagnetic_potential(m, rθ)
    -(utuϕ[1] + q * V[1])
end
function CircularOrbits.angmom(m::KerrNewmanMetric, rθ, utuϕ; q = 0.0)
    V = __KerrNewmanAD.electromagnetic_potential(m, rθ)
    (utuϕ[2] + q * V[4])
end

function CircularOrbits.Ω(
    m::KerrNewmanMetric,
    rθ;
    q = 0.0,
    μ = 1.0,
    contra_rotating = false,
    Ω_init = 1e-3,
)
    g, jacs = Gradus.metric_jacobian(m, rθ)
    # only want the derivatives w.r.t. r
    ∂rg = jacs[:, 1]

    # if q == 0.0 we have an analytic solution
    # since no maxwell tensor
    if q == 0.0
        return _Ω_analytic(∂rg, contra_rotating)
    end

    x = SVector(0, rθ[1], rθ[2], 0)
    F = Gradus.faraday_tensor(m, x)

    function f(ω)
        Δ = (ω^2 * ∂rg[4] + 2 * ω * ∂rg[5] + ∂rg[1])
        arg = -(ω^2 * g[4] + 2 * ω * g[5] + g[1]) / μ^2
        # analytic continuation
        ut = sign(arg) * sqrt(abs(arg))
        # combine with charge terms
        0.5 * Δ + (F[2, 4] * ω + F[2, 1]) * g[2] * q * ut
    end
    T = eltype(rθ)
    Ω₀ = (contra_rotating ? -T(Ω_init) : T(Ω_init))
    Roots.find_zero((f, w -> ForwardDiff.derivative(f, w)), Ω₀)
end

function find_isco_bounds(
    m::KerrNewmanMetric{T},
    max_upper_bound,
    step,
    ;
    kwargs...,
) where {T}
    # find upper bounds quick and dirty
    upper_bound = max_upper_bound
    for r = inner_radius(m):step:max_upper_bound
        # i don't want to have to try/catch here but
        # unsure how else to do this
        try
            CircularOrbits.energy(m, r; kwargs...)
        catch e
            if isa(e, Roots.ConvergenceFailed)
                break
            end
            rethrow(e)
        end
        upper_bound = r
    end
    # iterate in reverse with a negative step to find lower bound
    for r = upper_bound:(-step):1
        en = CircularOrbits.energy(m, r; kwargs...)
        if abs(en) > 1
            return r, upper_bound
        end
    end
    # for type stability
    return T(0), T(0)
end

export KerrNewmanMetric
