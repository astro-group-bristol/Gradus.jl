struct EmissivityProfileSetup{TryToOptimize,T,SamplerType,SpectrumType}
    λmax::T
    δmin::T
    δmax::T
    sampler::SamplerType
    spectrum::SpectrumType
    n_samples::Int
end

function EmissivityProfileSetup(
    T,
    spectrum;
    λmax = T(10000),
    δmin = T(0.01),
    δmax = T(179.99),
    sampler = nothing,
    n_samples = 1000,
    other_kwargs...,
)
    _sampler = if !isnothing(sampler)
        sampler
    else
        EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
    end
    setup = EmissivityProfileSetup{isnothing(sampler),T,typeof(_sampler),typeof(spectrum)}(
        λmax,
        δmin,
        δmax,
        _sampler,
        spectrum,
        n_samples,
    )
    other_kwargs, setup
end

"""
    source_to_disc_emissivity(
        m::AbstractStaticAxisSymmetric,
        spec::AbstractCoronalSpectrum,
        N,
        A,
        x,
        g,
        v_disc,
    )

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
    v_disc,
)
    # account for relativistic effects in area due to lorentz shift
    γ = lorentz_factor(m, x, v_disc)
    # divide by area to get number density
    I = coronal_spectrum(spec, g)
    N * I / (A * γ)
end

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
d = ThinDisc(Gradus.isco(m), 1000.0)
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
    m::AbstractMetric{T},
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel,
    spectrum = PowerLawSpectrum(2),
    ;
    kwargs...,
) where {T}
    solver_kwargs, setup = EmissivityProfileSetup(T, spectrum; kwargs...)
    emissivity_profile(setup, m, d, model; solver_kwargs...)
end

function emissivity_profile(
    setup::EmissivityProfileSetup,
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel,
    ;
    grid = GeometricGrid(),
    N = 100,
    kwargs...,
)
    RadialDiscProfile(
        tracecorona(
            m,
            d,
            model;
            sampler = setup.sampler,
            λmax = setup.λmax,
            n_samples = setup.n_samples,
        ),
        setup.spectrum;
        grid = grid,
        N = N,
    )
end

function _proper_area(m, x::SVector{4})
    gcomp = Gradus.metric_components(m, SVector(x[2], x[3]))
    det_g = √(gcomp[2] * gcomp[4])
    2π * det_g
end

function polar_angle_to_velfunc(m::AbstractMetric, x, v, δs; ϕ = zero(eltype(x)))
    function _polar_angle_velfunc(i)
        sky_angles_to_velocity(m, x, v, δs[i], ϕ)
    end
end

export emissivity_profile, tracecorona
