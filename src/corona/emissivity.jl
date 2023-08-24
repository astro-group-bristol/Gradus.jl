"""
    source_to_disc_emissivity(m, N, A, x, g, spec)

Compute the emissivity of a disc element with (proper) area `A` at coordinates `x` with metric
`m` and coronal spectrum `spec`. Since the emissivity is dependent on the incident flux, the photon (geodesic) count `N` must
be specified, along with the ratio of energies `g` (computed with [`energy_ratio`](@ref)) and the spectrum `spec`. 

The mathematical definition is
```math
\\varepsilon = \\frac{N}{A g^\\Gamma \\gamma}, 
```

where ``\\gamma`` is the Lorentz factor due to the velocity of the local disc frame. The velocity is currently
always considered to be the Keplerian velocity.

Wilkins & Fabian (2012) and Gonzalez et al. (2017).
"""
function source_to_disc_emissivity(
    m::AbstractStaticAxisSymmetric,
    spec::AbstractCoronalSpectrum,
    N,
    A,
    x,
    g,
)
    v = CircularOrbits.fourvelocity(m, SVector(x[2], x[3]))
    # account for relativistic effects in area due to lorentz shift
    γ = lorentz_factor(m, x, v)
    # divide by area to get number density
    I = coronal_spectrum(spec, g)
    N * I / (A * γ)
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
Jacobian of the spherical coordinates (with ``r`` fixes) yields a factor ``\\sin \\theta``
in order to maintain point density. It may be regarded as the PDF that samples ``\\theta`` uniformly.

Dauser et al. (2013)
"""
point_source_equatorial_disc_emissivity(spec::AbstractCoronalSpectrum, θ, g, A, γ) =
    sin(θ) * coronal_spectrum(spec, g) / (A * γ)

"""
    function emissivity_profile(
        m::AbstractMetric,
        d::AbstractAccretionGeometry,
        model::AbstractCoronaModel;
        kwargs...,
    end

Calculate the reflection emissivity profile of an accretion disc `d` around the spacetime `m`
for an illuminating coronal model `model`.

Returns a [`RadialDiscProfile`](@ref) via (Monte-Carlo or uniform) sampling
of the [`AbstractCoronaModel`](@ref) position and velocity distribution.

This function will attempt to automatically switch to use a better scheme to 
calculate the emissivity profiles if one is available. If not, the default 
algorithm is to estimate photon count ``N`` and calculate the emissivity with [`source_to_disc_emissivity`](@ref).

Common keyword arguments:
- `n_samples`: the maximum number of individual geodesics to sample on the emitter's sky.

Please consult the documentation of a specific model (e.g. [`LampPostModel`](@ref)) to
see algorithm specific keywords that may be passed.

All other keyword arguments are forwarded to [`tracegeodesics`](@ref).

## Example

```julia
m = KerrMetric()
d = GeometricThinDisc(Gradus.isco(m), 1000.0, π/2)
model = LampPostModel(h = 10.0)

profile = emissivity_profile(m, d, model; n_samples = 128)

# visualise as a function of disc radius
using Plots
plot(profile)
```

## Notes

The sampling is performed using an [`AbstractDirectionSampler`](@ref), 
which samples angles on the emitters sky along which a geodesic is traced. 
The effects of the spacetime and the observer's velocity are taken into account 
by using [`tetradframe`](@ref) and the corresponding coordinate transformation 
for local to global coordinates.

This function assumes axis symmetry, and therefore always interpolates the emissivity
as a function of the radial coordinate on the disc. If non-symmetric profiles are 
desired, consider using [`tracecorona`](@ref) with a profile constructor, e.g. [`VoronoiDiscProfile`](@ref).
"""
function emissivity_profile(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    spectrum = PowerLawSpectrum(2.0),
    kwargs...,
)
    emissivity_profile(m, d, model, spectrum; kwargs...)
end
function emissivity_profile(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel,
    spectrum::AbstractCoronalSpectrum;
    sampler = nothing,
    kwargs...,
)
    emissivity_profile(sampler, m, d, model, spectrum; kwargs...)
end
function emissivity_profile(
    sampler::AbstractDirectionSampler,
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel,
    spectrum::AbstractCoronalSpectrum;
    grid = GeometricGrid(),
    N = 100,
    kwargs...,
)
    RadialDiscProfile(
        tracecorona(m, d, model; sampler = sampler, kwargs...),
        spectrum;
        grid = grid,
        N = N,
    )
end
function emissivity_profile(
    ::Nothing,
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel,
    spectrum::AbstractCoronalSpectrum;
    kwargs...,
)
    @warn(
        "No sampler specified and no sophisticated algorithm for " *
        "$(Base.typename(typeof(model)).name) with $(Base.typename(typeof(d)).name) " *
        "implemented. Defaulting to `EvenSampler(BothHemispheres(), GoldenSpiralGenerator())` " *
        "as the sampler. Pass the `sampler` keyword to specify a different sampler."
    )
    emissivity_profile(
        EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
        m,
        d,
        model,
        spectrum;
        kwargs...,
    )
end

function _proper_area(m, x::SVector{4})
    gcomp = Gradus.metric_components(m, SVector(x[2], x[3]))
    det_g = √(gcomp[2] * gcomp[4])
    2π * det_g
end

function polar_angle_to_velfunc(m::AbstractMetric, x, v, δs)
    function _polar_angle_velfunc(i)
        sky_angles_to_velocity(m, x, v, δs[i], 0.0)
    end
end

function _point_source_symmetric_emissivity_profile(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel,
    spec::AbstractCoronalSpectrum;
    n_samples = 100,
    λ_max = 10_000.0,
    callback = domain_upper_hemisphere(),
    δmin = 0.01,
    δmax = 179.99,
    kwargs...,
)
    δs = deg2rad.(range(δmin, δmax, n_samples))
    # we assume a point source
    x, v = sample_position_velocity(m, model)
    velfunc = polar_angle_to_velfunc(m, x, v, δs)
    gps = tracegeodesics(
        m,
        x,
        velfunc,
        d,
        λ_max;
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
    J = sortperm(points, by = i -> _equatorial_project(i.x))
    points = points[J]
    δs = δs[J]

    r, ε = _point_source_emissivity(m, d, spec, v, δs, points)
    t = [i.x[1] for i in @views(points[1:end-1])]

    RadialDiscProfile(r, t, ε)
end

function _point_source_emissivity(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    spec::AbstractCoronalSpectrum,
    source_velocity,
    δs,
    points::AbstractVector{<:AbstractGeodesicPoint{T}},
) where {T}
    # function for obtaining keplerian velocities
    _disc_velocity = _keplerian_velocity_projector(m, d)
    # get radial coordinate
    r = [_equatorial_project(i.x) for i in points]

    # radial bin size
    Δr = diff(r)
    _points = @views(points[1:end-1])

    ε = map(enumerate(_points)) do (i, p)
        v_disc = _disc_velocity(p.x)
        gs = energy_ratio(m, p, source_velocity, v_disc)
        γ = lorentz_factor(m, p.x, v_disc)
        A = _proper_area(m, p.x) * Δr[i]

        point_source_equatorial_disc_emissivity(spec, δs[i], gs, A, γ)
    end

    r = r[1:end-1]
    r, ε
end

export emissivity_profile, tracecorona
