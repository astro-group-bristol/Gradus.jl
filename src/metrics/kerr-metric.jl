module __BoyerLindquistAD
using ..StaticArrays
using ..MuladdMacro

@muladd @fastmath begin
    Σ(r, a, θ) = r^2 + a^2 * cos(θ)^2
    Δ(r, R, a) = r^2 + a^2 - R * r

    # the way this function must be defined is a little complex
    # but helps with type-stability
    function metric_components(M, a, rθ)
        (r, θ) = rθ
        R = 2M
        sinθ2 = sin(θ)^2
        cosθ2 = (1 - sinθ2)
        # slightly faster, especially when considering AD evals
        Σ₀ = r^2 + a^2 * cosθ2

        tt = -(1 - (R * r) / Σ₀)
        rr = Σ₀ / Δ(r, R, a)
        θθ = Σ₀
        ϕϕ = sinθ2 * (r^2 + a^2 + (sinθ2 * R * r * a^2) / Σ₀)

        tϕ = (-R * r * a * sinθ2) / Σ₀
        @SVector [tt, rr, θθ, ϕϕ, tϕ]
    end
end

end # module

# new structure for our spacetime
"""
    struct KerrMetric{T} <: AbstractStaticAxisSymmetric{T}

The Kerr metric in Boyer-Lindquist coordinates, describing a black hole with mass ``M`` and
angular spin ``a``:

```math
\\begin{align*}
    \\text{d}s^2 =
    &- \\left( 1 - \\frac{2 M r}{\\Sigma} \\right)\\text{d}t^2
    - \\frac{2M r a \\sin^2(\\theta)}{\\Sigma} \\text{d}t \\text{d}\\phi
    \\\\
    &+ \\frac{\\Sigma}{\\Delta} \\text{d}r^2
    + \\Sigma \\text{d}\\theta^2
    + \\left(r^2 + a^2 + \\frac{2 M r a^2 \\sin^2(\\theta)}{\\Sigma} \\right) \\sin^2(\\theta) \\text{d}\\phi^2,
\\end{align*}
```
where
```math
\\Sigma = r^2 + a^2 \\cos^2 (\\theta),
\\quad \\text{and} \\quad
\\Delta = r^2 - 2Mr + a^2.
```

## Parameters
- `M = 1`: black hole mass.
- `a = 0`: black hole spin.
"""
@with_kw struct KerrMetric{T} <: AbstractStaticAxisSymmetric{T}
    @deftype T
    "Black hole mass."
    M = 1.0
    "Black hole spin."
    a = 0.0
end

# implementation
metric_components(m::KerrMetric, rθ) = __BoyerLindquistAD.metric_components(m.M, m.a, rθ)
inner_radius(m::KerrMetric) = m.M + √(m.M^2 - m.a^2)

# additional utilities
function convert_angles(a, r, θ, ϕ, θ_obs, ϕ_obs)
    ϕ̃ = ϕ - ϕ_obs
    R = sqrt(r^2 + a^2)
    o1 = r * R * sin(θ) * sin(θ_obs) * cos(ϕ̃) + R^2 * cos(θ) * cos(θ_obs)
    o2 = R * cos(θ) * sin(θ_obs) * cos(ϕ̃) - r * sin(θ) * cos(θ_obs)
    o3 = sin(θ_obs) * sin(ϕ̃) / sin(θ)
    sigma = r^2 + a^2 * cos(θ)^2
    -o1 / sigma, -o2 / sigma, o3 / R
end

# for disc profile models
function vector_to_local_sky(m::KerrMetric, u, θ, ϕ)
    convert_angles(m.a, u[2], u[3], u[4], θ, ϕ)
end

# special radii
isco(m::KerrMetric) = __BoyerLindquistFO.isco(m.M, m.a)

export KerrMetric
