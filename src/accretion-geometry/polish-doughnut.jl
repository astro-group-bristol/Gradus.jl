module __PolishDoughnut

import Gradus:
    CircularOrbits,
    AbstractMetric,
    __BoyerLindquistAD,
    KerrMetric,
    inverse_metric_components,
    metric_components
import ForwardDiff, Roots
import OrdinaryDiffEq: solve, Tsit5, DiscreteCallback, ODEProblem, terminate!
import StaticArrays: SVector

# modified Fuerst & Wu (2004, 2007)
function Ω(m::AbstractMetric, x, rₖ, n)
    rsinθ = x[1] * sin(x[2])
    xprime = SVector(rsinθ, π / 2)
    CircularOrbits.Ω(m, xprime) * (rₖ / (rsinθ))^n
end

function orbital_energy(m::AbstractMetric, x, rₖ, n)
    𝛺 = Ω(m, x, rₖ, n)
    ut_uϕ = CircularOrbits.ut_uϕ(𝛺, inverse_metric_components(m, x))
    CircularOrbits.energy(m, x, ut_uϕ)
    Ω₀ = 𝛺
    g = metric_components(m, x)
    -(g[1] + g[5] * Ω₀) / √(-g[1] - 2g[5] * Ω₀ - g[4] * Ω₀^2)
end

# Younsi et al. (2012): eq. (30)
Ψ₁(r, θ, M, a, invΩ, Σ) = M * ((Σ - 2r^2) / (Σ^2)) * (invΩ - a * sin(θ))^2 + r * sin(θ)^2
Ψ₁(m, x, 𝛺, Σ) = Ψ₁(x[1], x[2], m.M, m.a, 𝛺, Σ)

# Younsi et al. (2012): eq. (31)
Ψ₂(r, θ, M, a, invΩ, Σ, Δ) =
    sin(2θ) * ((M * r / (Σ^2)) * (a * invΩ - (r^2 + a^2))^2 + Δ / 2)
Ψ₂(m, x, 𝛺, Σ, Δ) = Ψ₂(x[1], x[2], m.M, m.a, 𝛺, Σ, Δ)

function isobar_differential(m::KerrMetric, x, inv𝛺)
    Σ = __BoyerLindquistAD.Σ(x[1], m.a, x[2])
    Δ = __BoyerLindquistAD.Δ(x[1], 2 * m.M, m.a)
    ψ₁ = Ψ₁(m, x, inv𝛺, Σ)
    ψ₂ = Ψ₂(m, x, inv𝛺, Σ, Δ)

    d = inv(√(Δ * ψ₁^2 + ψ₂^2) * √inv(Δ / Σ))
    dr = ψ₂ * d
    dθ = -ψ₁ * d

    SVector(dr, dθ)
end

function innermost_radius(m, rₖ, n; init_r = 5.0)
    # we want to find r for which dE/dr == 0
    f = r -> orbital_energy(m, SVector(r, π / 2), rₖ, n)
    df = r -> ForwardDiff.derivative(f, r)
    d2f = r -> ForwardDiff.derivative(df, r)
    Roots.find_zero((df, d2f), init_r)
end

function isobar(
    m::AbstractMetric,
    inner_radius,
    rₖ,
    n;
    λ_max = 40.0,
    dtmax = 5e-2,
    solver_opts...,
)
    terminate_when_negative =
        DiscreteCallback((u, t, integrator) -> u[1] * cos(u[2]) < 0, terminate!)

    u0 = SVector(inner_radius, π / 2)

    function drdθ(u, p, λ)
        𝛺 = Ω(m, u, rₖ, n)
        isobar_differential(m, u, inv(𝛺))
    end

    prob = ODEProblem{false}(drdθ, u0, (0.0, λ_max))
    sol = solve(
        prob,
        Tsit5();
        callback = terminate_when_negative,
        dtmax = dtmax,
        solver_opts...,
    )

    # unpack solution
    r = first.(sol.u)
    z = @. cos(last(sol.u)) * r

    # positivity
    I = @. z > 0
    r[I], z[I]
end

end # module

struct PolishDoughnut{T,F} <: AbstractThickAccretionDisc{T}
    inner_radius::T
    outer_radius::T
    rₖ::T
    n::T
    f::F
end

function PolishDoughnut(m::AbstractMetric; rₖ = 12.0, n = 0.21, init_r = 5.0, kwargs...)
    inner_radius = __PolishDoughnut.innermost_radius(m, rₖ, n; init_r = init_r)
    r, z = __PolishDoughnut.isobar(m, inner_radius, rₖ, n; kwargs...)
    # interpolate cross section
    f = DataInterpolations.LinearInterpolation(z, r)
    outer_radius = maximum(r)
    PolishDoughnut(inner_radius, outer_radius, rₖ, n, f)
end

function cross_section(d::PolishDoughnut, u4)
    let r = u4[2]
        if d.inner_radius ≤ r ≤ d.outer_radius
            d.f(r)
        else
            0.0
        end
    end
end

export PolishDoughnut
