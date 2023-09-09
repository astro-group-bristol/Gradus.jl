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

cross_section(d::ThickDisc, x4) = d.f(x4, d.params...)
cross_section(d::ThickDisc{T,F,Nothing}, x4) where {T,F} = d.f(x4)

xz_parameterize(d::AbstractAccretionDisc, ρ) =
    SVector(ρ, cross_section(d, SVector(0, ρ, π / 2, 0)))

function distance_to_disc(d::AbstractThickAccretionDisc, x4; gtol)
    height = cross_section(d, x4)
    if height <= 0
        return one(eltype(x4))
    end
    z = @inbounds x4[2] * cos(x4[3])
    abs(z) - height - (gtol * x4[2])
end

function _cartesian_tangent_vector(ρ, d::AbstractThickAccretionDisc)
    function _parameterization(r)
        xz_parameterize(d, r)
    end
    ∇f = ForwardDiff.derivative(_parameterization, ρ)
    v = SVector(∇f[1], 0, ∇f[2])
    v ./ √(v[1]^2 + v[2]^2 + v[3]^2)
end

function _cartesian_surface_normal(ρ, d::AbstractThickAccretionDisc)
    ∇f = _cartesian_tangent_vector(ρ, d)
    # rotate 90 degrees in about ϕ̂
    SVector(-∇f[3], ∇f[2], ∇f[1])
end

function _rotate_cartesian_about_z(n, ϕ)
    SVector(n[1] * cos(ϕ), n[1] * sin(ϕ), n[3])
end

function _cartesian_surface_normal(ρ, ϕ, d::AbstractThickAccretionDisc)
    _rotate_cartesian_about_z(_cartesian_surface_normal(ρ, d), ϕ)
end

export ThickDisc
