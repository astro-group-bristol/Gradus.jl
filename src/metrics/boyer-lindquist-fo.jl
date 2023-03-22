module __BoyerLindquistFO
using DocStringExtensions

"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.10):

```math
T = E \\left( r^2 + a^2 \\right) - L * a.
```
"""
@inline T(E, L, r, a) = E * (r^2 + a^2) - L * a


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.10), for the special case of a null-geodesic ``\\mu = 0``:

```math
V_r = T^2 - \\Delta \\left[ (L - a E)^2 + Q \\right]
```
"""
@inline Vr(E, L, M, Q, r, a) = T(E, L, r, a)^2 - Δ(M, r, a) * ((L - a * E)^2 + Q)


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.10), for the special case of a null-geodesic ``\\mu = 0``:

```math
V_\\theta =
    Q + \\cos^2 (\\theta) \\left[ a^2 E^2 - \\frac{L^2}{\\sin^2 (\\theta) } \\right].
```
"""
@inline Vθ(E, L, Q, a, θ) = Q + cos(θ)^2 * ((a * E)^2 - (L / sin(θ))^2)


"""
    $(TYPEDSIGNATURES)

The ``t`` compontent of the equation of motion for a photon around a black hole, multiplied
by ``\\Sigma``.

From Bardeen et al. (1972) eq. (2.9d):

```math
\\Sigma \\frac{\\text{d}t}{\\text{d}\\lambda} =
- a \\left( a E \\sin^2 \\theta - L \\right)
+ \\frac{\\left( r^2 + a^2 \\right) T}{\\Delta}.
```
"""
@inline function Σδt_δλ(E, L, M, r, a, θ)
    -a * (a * E * sin(θ)^2 - L) + (r^2 + a^2) * T(E, L, r, a) / Δ(M, r, a)
end


"""
    $(TYPEDSIGNATURES)

The ``r`` compontent of the equation of motion for a photon around a black hole, multiplied
by ``\\Sigma``.

Modified from Bardeen et al. (1972) eq. (2.9a):

```math
\\Sigma \\frac{\\text{d}r}{\\text{d}\\lambda} =
\\pm \\sqrt{\\lvert V_r \\rvert},
```

where, for implementation reason, the sign is always positive. Instead, the sign is applied
in [`δ`](@ref).
"""
@inline function Σδr_δλ(E, L, M, Q, r, a)
    V = Vr(E, L, M, Q, r, a)
    √abs(V)
end


"""
    $(TYPEDSIGNATURES)

The ``\\theta`` compontent of the equation of motion for a photon around a black hole,
multiplied by ``\\Sigma``.

Modified from Bardeen et al. (1972) eq. (2.9b):

```math
\\Sigma \\frac{\\text{d}\\theta}{\\text{d}\\lambda} =
\\pm \\sqrt{\\lvert V_\\theta \\rvert},
```

where, for implementation reason, the sign is always positive. Instead, the sign is applied
in [`δ`](@ref).
"""
@inline function Σδθ_δλ(E, L, Q, a, θ)
    V = Vθ(E, L, Q, a, θ)
    √abs(V)
end


"""
    $(TYPEDSIGNATURES)

The ``\\phi`` compontent of the equation of motion for a photon around a black hole,
multiplied by ``\\Sigma``.

From Bardeen et al. (1972) eq. (2.9c):

```math
\\Sigma \\frac{\\text{d}\\phi}{\\text{d}\\lambda} =
- \\frac{L}{\\sin^2 \\theta} - aE
+ \\frac{aT}{\\Delta}.
```
"""
@inline function Σδϕ_δλ(E, L, M, r, a, θ)
    (L / sin(θ)^2) - (a * E) + a * T(E, L, r, a) / Δ(M, r, a)
end


"""
    $(TYPEDSIGNATURES)

From Fanton et al. (1997), eq. (76):

```math
S = 1 + \\alpha \\omega \\sin \\theta.
```
"""
S(θ, α, ω) = 1 + α * ω * sin(θ)


"""
    $(TYPEDSIGNATURES)

Calculates and returns the observer's angles ``\\sin \\Theta`` and ``\\sin \\Phi``, where
the parameters in the function signature correspond to the Bardeen et al. (1972) paper.

From Fanton et al. (1997), eq. (74) and eq. (75):

```math
\\begin{align*}
\\sin \\Theta &=
\\frac{\\alpha \\Sigma}{\\sqrt{A}}
\\left\\{
        \\beta^2 + \\left( \\alpha + a \\sin \\theta \\right)^2
        + \\frac{
            A S^2 - \\left( r^2 + a^2 + a \\alpha \\sin \\theta \\right)
        }{\\Delta}
\\right\\}, \\\\
\\sin \\Phi &=
-\\frac{\\alpha \\Sigma \\sqrt{\\Delta}}{A S \\sin\\Theta}.
\\end{align*}
```
"""
function sinΦsinΨ(Σ₀, sinθ, A₀, Δ₀, S₀, r, a, α, β)
    # calc 1
    sinΦ =
        ((α * Σ₀) / √A₀) /
        sqrt(β^2 + (α + a * sinθ)^2 + (A₀ * S₀^2 - (r^2 + a^2 + a * α * sinθ)^2) / Δ₀)

    # calc 2
    sinΨ = -(α * Σ₀ * √Δ₀) / (S₀ * A₀ * sinΦ)

    # return
    (sinΦ, sinΨ)
end

function sinΦsinΨ(M, r, a, θ, α, β)
    # value cache
    Σ₀ = Σ(r, a, θ)
    sinθ = sin(θ)
    A₀ = A(M, r, a, θ)
    Δ₀ = Δ(M, r, a)
    ω = 2.0 * a * r / A₀
    S₀ = S(θ, α, ω)

    sinΦsinΨ(Σ₀, sinθ, A₀, Δ₀, S₀, r, a, α, β)
end


"""
    $(TYPEDSIGNATURES)

Calculates conserved quantities
- angular momentum ``L``
- Carter parameter ``Q``

for a photon described with position described by `x` in a Kerr spacetime given by
`p`.

From Fanton et al. (1997), eq. (69):

```math
L = \\frac{\\Upsilon_1}{\\Upsilon_2},
```

where

```math
\\begin{align*}
\\Upsilon_1 &= \\sin \\theta \\sin \\Phi \\sin \\Theta, \\\\
\\Upsilon_2 &= \\frac{\\Sigma \\sqrt{\\Delta}}{A} + \\omega \\Upsilon_1,
\\omega = \\frac{2 a r}{A},
\\end{align*}
```

taken from eq. (72) and eq. (73).

From Fanton et al. (1997), eq. (70):

```math
Q =
\\frac{P^2}{\\Delta}
- \\left( \\lambda - a \\right)^2
- \\frac{\\Sigma^2}{A} \\left( \\frac{\\cos \\Phi}{\\Upsilon_2} \\right)^2,
```

and

```math
P = \\left( r^2 + a^2 - a L \\right),
```

taken from eq. (71).
"""
function LQ(M, r, a, θ, sinΦ, sinΨ)
    sinθ = sin(θ)
    Σ₀ = Σ(r, a, θ)
    A₀ = A(M, r, a, θ)
    Δ₀ = Δ(M, r, a)
    ω = 2.0 * a * r / A₀

    Υ₁ = sinθ * sinΦ * sinΨ
    Υ₂ = (Σ₀ * √Δ₀ / A₀) + ω * Υ₁

    Υ₁ = sinθ * sinΦ * sinΨ
    Υ₂ = (Σ₀ * √Δ₀ / A₀) + ω * Υ₁

    L = Υ₁ / Υ₂
    P = (r^2 + a^2 - (a * L))
    Q = (P^2 / Δ₀) - (L - a)^2 - ((Σ₀^2 / A₀) * (cos(asin(sinΨ)) / Υ₂)^2)
    L, Q
end


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.3):

```math
\\Sigma = r^2 + a^2 \\cos^2( \\theta ).
```
"""
Σ(r, a, θ) = r^2 + (a * cos(θ))^2


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.3):

```math
\\Delta = r^2 - 2 M r + a^2.
```
"""
Δ(M, r, a) = r^2 - 2 * M * r + a^2


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.3):

```math
A = (r^2 + a^2)^2 - a^2 \\Delta \\sin^2 ( \\theta ).
```
"""
A(M, r, a, θ) = (r^2 + a^2)^2 - a^2 * Δ(M, r, a) * sin(θ)^2


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.21):

```math
Z_1 =
1 + \\sqrt[3]{1 - \\frac{a^2}{M^2}}
\\left[
    \\sqrt[3]{1 + \\frac{a}{M}} + \\sqrt[3]{1 - \\frac{a}{M}}
\\right].
```
"""
Z₁(M, a) = 1 + ∛(1 - (a / M)^2) * (∛(1 + (a / M)) + ∛(1 - (a / M)))


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.21):

```math
Z_2 = \\sqrt{\\frac{3a^2}{M^2} + Z_1^2}.
```
"""
Z₂(M, a) = √(3(a / M)^2 + Z₁(M, a)^2)


function four_velocity(u, E, M, a, p)
    let r = u[2], θ = u[3], L = p.L, Q = p.Q
        Σ₀ = Σ(r, a, θ)
        (
            Σδt_δλ(E, L, M, r, a, θ) / Σ₀,
            p.r * Σδr_δλ(E, L, M, Q, r, a) / Σ₀,
            p.θ * Σδθ_δλ(E, L, Q, a, θ) / Σ₀,
            Σδϕ_δλ(E, L, M, r, a, θ) / Σ₀,
        )
    end
end

"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.21):

```math
r_\\text{ms} = M \\left\\{ 3 + Z_2 \\pm \\sqrt{(3 - Z_1)(3 + Z_1 + 2 Z_2)} \\right\\}.
```

The choice of ``\\pm`` is chosen by the sign of ``a``.
"""
isco(M, a, ±) = M * (3 + Z₂(M, a) ± √((3 - Z₁(M, a)) * (3 + Z₁(M, a) + 2 * Z₂(M, a))))
isco(M, a) = a > 0.0 ? isco(M, a, -) : isco(M, a, +)
end # module

"""
A first-order implementation of [`KerrMetric`](@ref).
$(FIELDS)
"""
@with_kw struct KerrSpacetimeFirstOrder{T} <: AbstractFirstOrderMetricParams{T}
    @deftype T
    "Black hole mass."
    M = 1.0
    "Black hole spin."
    a = 0.0
    "Geodesic energy (a consant of motion)."
    E = 1.0
end

GradusBase.inner_radius(m::KerrSpacetimeFirstOrder{T}) where {T} = m.M + √(m.M^2 - m.a^2)
GradusBase.constrain(::KerrSpacetimeFirstOrder{T}, u, v; μ = 0.0) where {T} = v[1]

four_velocity(u, m::KerrSpacetimeFirstOrder{T}, p) where {T} =
    __BoyerLindquistFO.four_velocity(u, m.E, m.M, m.a, p)
calc_lq(m::KerrSpacetimeFirstOrder{T}, pos, vel) where {T} =
    __BoyerLindquistFO.LQ(m.M, pos[2], m.a, pos[3], vel[3], vel[4])

function Vr(m::KerrSpacetimeFirstOrder{T}, u, p) where {T}
    L = p.L
    Q = p.Q
    __BoyerLindquistFO.Vr(m.E, L, m.M, Q, u[2], m.a)
end
function Vθ(m::KerrSpacetimeFirstOrder{T}, u, p) where {T}
    L = p.L
    Q = p.Q
    __BoyerLindquistFO.Vθ(m.E, L, Q, m.a, u[3])
end

function impact_parameters_to_vel(m::KerrSpacetimeFirstOrder{T}, u, α, β) where {T}
    sinΦ, sinΨ = __BoyerLindquistFO.sinΦsinΨ(m.M, u[2], m.a, u[3], α, β)
    (β < 0.0 ? 1.0 : -1.0, sinΦ, sinΨ)
end

# special radii
isco(m::KerrSpacetimeFirstOrder{T}) where {T} = __BoyerLindquistFO.isco(m.M, m.a)

export KerrSpacetimeFirstOrder
