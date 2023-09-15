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

ThickDisc(::Number, ::Number) =
    error("Invalid constructor (you probably meant ThinDisc, not ThickDisc).")

function ThickDisc(cross_section::F) where {F}
    # todo: float bit generic???
    ThickDisc{Float64,F,Nothing}(cross_section, nothing)
end

function ThickDisc(cross_section::F, parameters::P) where {F,P}
    # todo: float bit generic???
    ThickDisc{Float64,F,P}(cross_section, parameters)
end

cross_section(d::ThickDisc, ρ) = d.f(ρ, d.params...)
cross_section(d::ThickDisc{T,F,Nothing}, ρ) where {T,F} = d.f(ρ)

xz_parameterize(d::AbstractAccretionDisc, ρ) = SVector(ρ, cross_section(d, ρ))

function distance_to_disc(d::AbstractThickAccretionDisc, x4; gtol)
    height = cross_section(d, _equatorial_project(x4))
    if height <= 0
        return one(eltype(x4))
    end
    _spinaxis_project(x4) - height
end

function _cartesian_tangent_vector(d::AbstractThickAccretionDisc, ρ)
    function _parameterization(ρ̄)
        xz_parameterize(d, ρ̄)
    end
    ∇f = ForwardDiff.derivative(_parameterization, ρ)
    v = SVector(∇f[1], 0, ∇f[2])
    v ./ √(v[1]^2 + v[2]^2 + v[3]^2)
end

function _cartesian_surface_normal(d::AbstractThickAccretionDisc, ρ)
    ∇f = _cartesian_tangent_vector(d, ρ)
    # rotate 90 degrees in about ϕ̂
    SVector(-∇f[3], ∇f[2], ∇f[1])
end

_cartesian_surface_normal(d::AbstractThickAccretionDisc, ρ, ϕ) =
    _rotate_about_spinaxis(_cartesian_surface_normal(d, ρ), ϕ)

export ThickDisc
