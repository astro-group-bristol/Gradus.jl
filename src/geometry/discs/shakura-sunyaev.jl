"""
    ShakuraSunyaev{T} <: AbstractThickAccretionDisc{T}
    ShakuraSunyaev(
        m::AbstractMetric;
        eddington_ratio = 0.3,
        η = nothing,
        contra_rotating = false,
    )

The classic Shakura & Sunyaev (1973) accretion disc model, with height given by ``2H``, where

```math
H = \\frac{3}{2} \\frac{1}{\\eta} \\left( \\frac{\\dot{M}}{\\dot{M}_\\text{Edd}} \\right) \\left( 1 - \\sqrt{\\frac{r_\\text{isco}}{\\rho}} \\right)
```

Here ``\\eta`` is the radiative efficiency, which, if unspecified, is determined by the circular orbit energy at the ISCO:

```math
\\eta = 1 - E_\\text{isco}
```
"""
struct ShakuraSunyaev{T} <: AbstractThickAccretionDisc{T}
    Ṁ_Ṁedd::T
    inv_η::T
    inner_radius::T
end

function cross_section(d::ShakuraSunyaev, ρ)
    if ρ < d.inner_radius
        return -zero(typeof(ρ))
    end
    3 * d.inv_η * d.Ṁ_Ṁedd * (1 - sqrt(d.inner_radius / ρ))
end

function ShakuraSunyaev(
    m::AbstractMetric{T};
    eddington_ratio = 0.3,
    η = nothing,
    contra_rotating = false,
) where {T}
    r_isco = isco(m)
    radiative_efficiency = if isnothing(η)
        1 - CircularOrbits.energy(
            m,
            SVector{2}(r_isco, π / 2);
            contra_rotating = contra_rotating,
        )
    else
        η
    end
    ShakuraSunyaev(T(eddington_ratio), inv(radiative_efficiency), r_isco)
end

export ShakuraSunyaev
