function build_collision_callback(
    g::AbstractAccretionDisc{T};
    gtol,
    interp_points = 8,
) where {T}
    ContinuousCallback(
        (u, λ, integrator) -> distance_to_disc(g, u; gtol = gtol),
        i -> terminate!(i, :Intersected);
        interp_points = interp_points,
        save_positions = (true, false),
    )
end

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

@fastmath function distance_to_disc(d::GeometricThinDisc{T}, u4; gtol) where {T}
    p = @inbounds let r = u4[2], θ = u4[3], ϕ = u4[4]
        if r < d.inner_radius || r > d.outer_radius
            return 1.0
        end
        sinθ = sin(θ)
        @SVector [r * sinθ * cos(ϕ), r * sinθ * sin(ϕ), r * cos(θ)]
    end
    n = @SVector [T(0.0), cos(d.inclination), sin(d.inclination)]
    # project u into normal vector n
    k = p ⋅ n
    abs(k) - (gtol * u4[2])
end

abstract type AbstractThickAccretionDisc{T} <: AbstractAccretionDisc{T} end
cross_section(d::AbstractThickAccretionDisc, u4) =
    error("Not implemented for $(typeof(d)).")

@with_kw struct ThickDisc{T,F,P} <: AbstractThickAccretionDisc{T}
    f::F
    params::P
end

function ThickDisc(cross_section::F, parameters::P; T = Float64) where {F,P}
    ThickDisc{T,F,P}(cross_section, parameters)
end

cross_section(d::ThickDisc, u4) = d.f(u4, d.params...)

function distance_to_disc(d::AbstractThickAccretionDisc, u4; gtol)
    height = cross_section(d, u4)
    if height <= 0.0
        return 1.0
    end
    z = @inbounds u4[2] * cos(u4[3])
    abs(z) - height - (gtol * u4[2])
end

# common thick disc models

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
    m::AbstractMetricParams{T};
    eddington_ratio = 0.3,
    η = nothing,
    contra_rotating = false,
) where {T}
    r_isco = isco(m)
    radiative_efficieny = if isnothing(η)
        1 - CircularOrbits.energy(
            m,
            SVector{2}(r_isco, π / 2);
            contra_rotating = contra_rotating,
        )
    else
        η
    end
    ShakuraSunyaev(T(eddington_ratio), 1.0, radiative_efficieny, r_isco)
end

export AbstractThickAccretionDisc, ThickDisc, ShakuraSunyaev, GeometricThinDisc
