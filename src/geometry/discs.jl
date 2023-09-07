"""
    distance_to_disc(d::AbstractAccretionGeometry, u; kwargs...)

Calculate distance to the closest element of the disc. The distance need not be metric or Pythagorean, but rather should be positive
when the four vector `u` is distant, zero when `u` is on the surface, and negative when `u` is within the disc geometry.

Must return a floating point number.
"""
function distance_to_disc(d::AbstractAccretionGeometry, u; kwargs...)
    if is_finite_disc(d)
        error("Not implemented for $(typeof(d)).")
    else
        1
    end
end


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
    # project u into equatorial plane
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

include("discs/datum-plane.jl")
include("discs/polish-doughnut.jl")
include("discs/shakura-sunyaev.jl")
include("discs/thick-disc.jl")
include("discs/thin-disc.jl")

export AbstractThickAccretionDisc, cross_section, PrecessingDisc, EllipticalDisc
