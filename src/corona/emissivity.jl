"""
    source_to_disc_emissivity(m, 𝒩, A, x, g)

Compute the emissivity (in arbitrary units) in the disc element with proper area `A`, photon
count `𝒩`, central position `x`, and redshift `g`. Evaluates

```math
\\varepsilon = \\frac{\\mathscr{N}}{A g^2}.
```
"""
function source_to_disc_emissivity(m::AbstractStaticAxisSymmetric, 𝒩, A, x, g; Γ = 2)
    v = CircularOrbits.fourvelocity(m, SVector(x[2], x[3]))
    # account for relativistic effects in area due to lorentz shift
    γ = lorentz_factor(m, x, v)
    # divide by area to get number density
    𝒩 / (g^Γ * A * γ)
end

function source_to_disc_emissivity(δ, g, A, γ; Γ = 2)
    sin(δ) / (g^Γ * A * γ)
end

"""
function emissivity_profile(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    N = 100,
    kwargs...,
end

Calculates a [`RadialDiscProfile`](@ref).

This function assumes axis symmetry, and therefore always interpolates the emissivity
as a function of the radial coordinate on the disc. If non-symmetric profiles are 
desired, consider using [`VoronoiDiscProfile`](@ref).
"""
emissivity_profile(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    sampler = nothing,
    kwargs...,
) = emissivity_profile(sampler, m, d, model; kwargs...)
function emissivity_profile(
    sampler::AbstractDirectionSampler,
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    grid = GeometricGrid(),
    N = 100,
    kwargs...,
)
    RadialDiscProfile(
        tracecorona(m, d, model; sampler = sampler, kwargs...);
        grid = grid,
        N = N,
    )
end
function emissivity_profile(
    ::Nothing,
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    kwargs...,
)
    error(
        """Not yet implemented for $(Base.typename(typeof(model)).name) with $(Base.typename(typeof(d)).name). 
        This dispatch is reserved for more sophisticated (rather that just pure Monte-Carlo sampling) 
        of the emissivity profile, and may not be applicable for all coronal models. Pass the `sampler` kwarg
        with an `AbstractDirectionSampler` to use the general strategy.
        """,
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
    model::AbstractCoronaModel;
    n_samples = 100,
    λ_max = 10_000.0,
    callback = domain_upper_hemisphere(),
    δmin = 0.1,
    δmax = 179.9,
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
    J = sortperm(points, by = i -> i.x[2])
    points = points[J]
    δs = δs[J]

    r, ε = _point_source_emissivity(m, d, v, δs, points)
    t = [i.x[1] for i in @views(points[1:end-1])]

    RadialDiscProfile(r, t, ε)
end

function _point_source_emissivity(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    source_velocity,
    δs,
    points,
)
    _points = @views points[1:end-1]
    # get the time and radial position
    r = [i.x[2] for i in points]
    # create a view
    gs = Gradus.energy_ratio.(m, _points, (source_velocity,))
    A = [_proper_area(m, i.x) for i in _points]
    γ = [Gradus.lorentz_factor(m, d, i.x) for i in _points]
    Δr = diff(r)
    @. A = A * Δr
    r = r[1:end-1]
    ε = source_to_disc_emissivity.(@views(δs[1:end-1]), gs, A, γ)
    r, ε
end

export emissivity_profile, tracecorona
