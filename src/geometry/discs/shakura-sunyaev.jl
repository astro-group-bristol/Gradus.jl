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
    Ṁ::T
    Ṁedd::T
    η::T
    risco::T
end

@fastmath function cross_section(d::ShakuraSunyaev, ρ)
    if ρ < d.risco
        return -one(typeof(ρ))
    end
    H = (3 / 2) * inv(d.η) * (d.Ṁ / d.Ṁedd) * (1 - sqrt(d.risco / ρ))
    2H
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
    ShakuraSunyaev(T(eddington_ratio), 1.0, radiative_efficiency, r_isco)
end

export ShakuraSunyaev
