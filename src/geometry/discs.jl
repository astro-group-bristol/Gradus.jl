"""
    distance_to_disc(d::AbstractAccretionGeometry, u; kwargs...)

Calculate distance to the closest element of the disc. The distance need not be metric or Pythagorean, but rather should be positive
when the four vector `u` is distant, zero when `u` is on the surface, and negative when `u` is within the disc geometry.

Must return a floating point number.
"""
distance_to_disc(d::AbstractAccretionGeometry, u; kwargs...) =
    error("Not implemented for $(typeof(d)).")

"""
    cross_section(d::AbstractThickAccretionDisc, u)

Return the height cross-section of a thick accretion disc at the (projected) coordinates of `u`. This function also incorporates bounds checking, and should
return a negative value if the disc is not defined at `u`.

## Example

For a top hat disc profile with constant height between two radii

```julia
struct TopHatDisc{T} <: AbstractThickAccretionDisc{T}
    inner_r::T
    outer_r::T
end

function Gradus.cross_section(d::TopHatDisc, u)
    # project u into equitorial plane
    r = u[2] * sin(u[3])
    if (r < d.inner_r) || (r > d.outer_r)
        return -1.0
    else
        return 1.0
    end
end
```
"""
cross_section(d::AbstractThickAccretionDisc, x4) =
    error("Not implemented for $(typeof(d)).")
r_cross_section(d::AbstractThickAccretionDisc, r::Number) =
    cross_section(d, SVector(0.0, r, π / 2, 0.0))

"""
    struct GeometricThinDisc{T} <: AbstractAccretionDisc{T}
    GeometricThinDisc(inner_radius::T, outer_radius::T, inclination::T)

$(FIELDS)

Simple geometrically thin accretion disc spanning from `inner_radius` to `outer_radius` in
gravitational units. Inclination of the disc is relative to spin axis, with ``90^\\circ``
being perpendicular to the spin axis.
"""
@with_kw struct GeometricThinDisc{T} <: AbstractAccretionDisc{T}
    inner_radius::T
    outer_radius::T
    inclination::T
end
optical_property(::Type{<:GeometricThinDisc}) = OpticallyThin()

@fastmath function distance_to_disc(d::GeometricThinDisc{T}, x4; gtol) where {T}
    p = @inbounds let r = x4[2], θ = x4[3], ϕ = x4[4]
        if r < d.inner_radius || r > d.outer_radius
            return 1.0
        end
        sinθ = sin(θ)
        @SVector [r * sinθ * cos(ϕ), r * sinθ * sin(ϕ), r * cos(θ)]
    end
    n = @SVector [T(0.0), cos(d.inclination), sin(d.inclination)]
    # project u into normal vector n
    k = p ⋅ n
    abs(k) - (gtol * x4[2])
end

"""
    ThickDisc{T,F,P} <: AbstractThickAccretionDisc{T}
    ThickDisc(f, params=nothing; T = Float64)

A standard wrapper for creating custom disc profiles from height cross-section function `f`. This function
is given the disc parameters as unpacked arguments:

```julia
d.f(u, d.params...)
```

If no parameters are specified, none will be passed to `f`.

## Example

Specifying a toroidal disc centered on `r=10` with radius 1:

```julia
d = ThickDisc() do u
    r = u[2]
    if r < 9.0 || r > 11.0
        return -1.0
    else
        x = r - 10.0
        sqrt(1-x^2)
    end
end
```

"""
struct ThickDisc{T,F,P} <: AbstractThickAccretionDisc{T}
    f::F
    params::P
end

function ThickDisc(cross_section::F) where {F}
    # todo: float bit generic???
    ThickDisc{Float64,F,Nothing}(cross_section, nothing)
end

function ThickDisc(cross_section::F, parameters::P) where {F,P}
    # todo: float bit generic???
    ThickDisc{Float64,F,P}(cross_section, parameters)
end

cross_section(d::ThickDisc, x4) = d.f(x4, d.params...)
cross_section(d::ThickDisc{T,F,Nothing}, x4) where {T,F} = d.f(x4)

function distance_to_disc(d::AbstractThickAccretionDisc, x4; gtol)
    height = cross_section(d, x4)
    if height <= 0.0
        return 1.0
    end
    z = @inbounds x4[2] * cos(x4[3])
    abs(z) - height - (gtol * x4[2])
end

# common thick disc models

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

@fastmath function cross_section(d::ShakuraSunyaev, u)
    @inbounds let r = u[2], θ = u[3]
        if r < d.risco
            return -1.0
        end
        ρ = r * sin(θ)
        H = (3 / 2) * inv(d.η) * (d.Ṁ / d.Ṁedd) * (1 - sqrt(d.risco / ρ))
        2H
    end
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

struct EllipticalDisc{T} <: AbstractAccretionDisc{T}
    inner_radius::T
    semi_major::T
    semi_minor::T
end

function distance_to_disc(d::EllipticalDisc, x4; gtol)
    if d.semi_major < x4[2] || x4[2] < d.inner_radius
        return 1.0
    end
    # equation of ellipse
    y = √((1 - (x4[2] / d.semi_major)^2) * d.semi_minor^2)
    # check height less than y with tolerance
    h = abs(x4[2] * cos(x4[3]))
    h - y - (gtol * x4[2])
end

struct PrecessingDisc{T,D<:AbstractAccretionDisc{T}} <: AbstractAccretionDisc{T}
    disc::D
    β::T
    γ::T
    R::SMatrix{3,3,T,9}
end

function PrecessingDisc(disc, β, γ)
    Rx = SMatrix{3,3}(1, 0, 0, 0, cos(-β), -sin(-β), 0, sin(-β), cos(-β))
    R = Rx
    PrecessingDisc(disc, β, γ, R)
end

function distance_to_disc(d::PrecessingDisc, x4; gtol)
    x = let θ = x4[3]
        ϕ = x4[4] - d.γ
        d.R * SVector(sin(θ) * sin(ϕ), sin(θ) * cos(ϕ), cos(θ))
    end
    x4prime = SVector(x4[1], x4[2], atan(√(x[1]^2 + x[2]^2), x[3]), atan(x[2], x[1]))
    distance_to_disc(d.disc, x4prime; gtol)
end

export AbstractThickAccretionDisc,
    ThickDisc,
    ShakuraSunyaev,
    GeometricThinDisc,
    cross_section,
    PrecessingDisc,
    EllipticalDisc
