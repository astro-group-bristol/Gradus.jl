@with_kw struct LampPostModel{T} <: AbstractCoronaModel{T}
    @deftype T
    h = 5.0
    θ = 0.01
    ϕ = 0.0
end

function sample_position_velocity(m::AbstractMetric, model::LampPostModel{T}) where {T}
    x = SVector{4,T}(0, model.h, model.θ, model.ϕ)
    gcomp = metric_components(m, SVector(x[2], x[3]))
    v = inv(√(-gcomp[1])) * SVector{4,T}(1, 0, 0, 0)
    x, v
end

"""
struct BeamedPointSource{T} <: Gradus.AbstractCoronaModel{T}
    r::T
    β::T
end

Point source corona moving away from the black hole with specified starting height above the disc `r` and source speed `β`.
Point sources are the most plausible example of source that would support beaming (ref: Gonzalez et al 2017).

"""
struct BeamedPointSource{T} <: Gradus.AbstractCoronaModel{T}
    r::T
    β::T
end

# these are specific to BeamedPointSource, so we'll scope them in a module
# so that they don't pollute the namespace of Gradus
module __BeamedPointSource
import ..Gradus: AbstractMetric, metric_components, SVector

drdt(g, β) = β * √(-g[1] / g[2])
drdt(m::AbstractMetric, x, β) = drdt(metric_components(m, SVector(x[2], x[3])), β)
end # module

function sample_position_velocity(m::AbstractMetric, model::BeamedPointSource)
    x = SVector{4}(0, model.r, 1e-4, 0)
    g = metric_components(m, SVector(x[2], x[3]))
    v̄ = SVector(1, __BeamedPointSource.drdt(g, model.β), 0, 0)
    v = constrain_normalize(m, x, v̄; μ = 1)
    x, v
end

"""
    point_source_equatorial_disc_emissivity(θ, g, A, γ, spec)

Calculate the emissivity of a point illuminating source on the spin axis for an annulus of the
equatorial accretion disc with (proper) area `A`. The precise formulation follows from Dauser et al. (2013),
with the emissivity calculated as
```math
\\varepsilon = \\frac{\\sin \\theta}{A g^\\Gamma \\gamma}
```
where ``\\gamma`` is the Lorentz factor due to the velocity of the local disc frame. 
The ratio of energies is `g` (computed with [`energy_ratio`](@ref)), with `spec` being the abstract
coronal spectrum and  ``\\theta`` is the angle from the spin axis in the emitters from at which the geodesic was directed. 
`coronal_spectrum` function is used to calculate the spectrum of the corona by taking `g` to the power of `Γ`, allowing
further modification of spectrum if needed, based on the value of the photon index.

The ``\\sin \\theta`` term appears to extend the result to three dimensions, since the
Jacobian of the spherical coordinates (with ``r`` fixed) yields a factor ``\\sin \\theta``
in order to maintain point density. It may be regarded as the PDF that samples ``\\theta`` uniformly.

Dauser et al. (2013)
"""
function point_source_equatorial_disc_emissivity(spec::AbstractCoronalSpectrum, θ, g, A, γ)
    sin(θ) * coronal_spectrum(spec, g) / (A * γ)
end


function _point_source_symmetric_emissivity_profile(
    setup::EmissivityProfileSetup,
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    callback = domain_upper_hemisphere(),
    kwargs...,
)
    δs = deg2rad.(range(setup.δmin, setup.δmax, setup.n_samples))
    # we assume a point source
    x, v = sample_position_velocity(m, model)
    velfunc = polar_angle_to_velfunc(m, x, v, δs)
    gps = tracegeodesics(
        m,
        x,
        velfunc,
        d,
        setup.λmax;
        save_on = false,
        ensemble = EnsembleEndpointThreads(),
        callback = callback,
        trajectories = length(δs),
        kwargs...,
    )

    # filter only those that intersected, and sort radially
    I = [i.status == StatusCodes.IntersectedWithGeometry for i in gps]
    points = gps[I]
    δs = δs[I]
    rs = [_equatorial_project(i.x) for i in points]
    J = sortperm(rs)
    rs_sorted = rs[J]
    points = points[J]
    δs = δs[J]

    r, ε = _point_source_emissivity(m, d, setup.spectrum, v, rs_sorted, δs, points)
    t = [i.x[1] for i in @views(points[1:end-1])]

    RadialDiscProfile(r, t, ε)
end

function _point_source_emissivity(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    spec::AbstractCoronalSpectrum,
    source_velocity,
    r,
    δs,
    points::AbstractVector{<:AbstractGeodesicPoint{T}},
) where {T}
    # function for obtaining keplerian velocities
    _disc_velocity = _keplerian_velocity_projector(m, d)

    # radial bin size
    _points = @views(points[1:end-1])

    ε = map(enumerate(_points)) do (i, p)
        v_disc = _disc_velocity(p.x)
        gs = energy_ratio(m, p, source_velocity, v_disc)
        γ = lorentz_factor(m, p.x, v_disc)
        Δr = abs(r[i+1] - r[i])
        A = _proper_area(m, p.x) * Δr

        point_source_equatorial_disc_emissivity(spec, δs[i], gs, A, γ)
    end

    r = r[1:end-1]
    r, ε
end

function emissivity_profile(
    setup::EmissivityProfileSetup{true},
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::Union{LampPostModel,BeamedPointSource};
    kwargs...,
)
    _point_source_symmetric_emissivity_profile(setup, m, d, model; kwargs...)
end

export LampPostModel, BeamedPointSource
