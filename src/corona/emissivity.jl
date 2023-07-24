"""
    source_to_disc_emissivity(m, 𝒩, A, x, g)

Compute the emissivity (in arbitrary units) in the disc element with proper area `A`, photon
count `𝒩`, central position `x`, and redshift `g`. Evaluates

```math
\\varepsilon = \\frac{\\mathscr{N}}{A g^2}.
```
"""
function source_to_disc_emissivity(m::AbstractStaticAxisSymmetric, 𝒩, A, x, g; Γ = 2)
    v = CircularOrbits.fourvelocity(m, x)
    # account for relativistic effects in area due to lorentz shift
    γ = lorentz_factor(m, SVector(0, x[1], x[2], 0), v)
    # divide by area to get number density
    𝒩 / (g^Γ * A * γ)
end

function source_to_disc_emissivity(δ, g, A, γ; Γ = 2)
    sin(δ) / (g^Γ * A * γ)
end

struct CoronaGeodesics{T,M,G,C,P,V}
    trace::T
    metric::M
    geometry::G
    model::C
    geodesic_points::P
    source_velocity::V
end

function tracecorona(
    m::AbstractMetric,
    g::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    λ_max = 10_000,
    n_samples = 1024,
    sampler = EvenSampler(domain = BothHemispheres(), generator = RandomGenerator()),
    trace = TraceGeodesic(),
    callback = domain_upper_hemisphere(),
    kwargs...,
)
    xs, vs, source_vels = sample_position_direction_velocity(m, model, sampler, n_samples)
    gps = tracegeodesics(
        m,
        xs,
        vs,
        g,
        λ_max;
        trace = trace,
        save_on = false,
        ensemble = EnsembleEndpointThreads(),
        callback = callback,
        kwargs...,
    )
    mask = [i.status == StatusCodes.IntersectedWithGeometry for i in gps]
    CoronaGeodesics(trace, m, g, model, gps[mask], source_vels[mask])
end

function RadialDiscProfile(ce::CoronaGeodesics; kwargs...)
    J = sortperm(ce.geodesic_points; by = i -> i.x[2])
    @views RadialDiscProfile(
        ce.metric,
        ce.model,
        ce.geodesic_points[J],
        ce.source_velocity[J],
    )
end

function RadialDiscProfile(ce::CoronaGeodesics{<:TraceRadiativeTransfer}; kwargs...)
    J = sortperm(ce.geodesic_points; by = i -> i.x[2])
    @views RadialDiscProfile(
        ce.metric,
        ce.model,
        ce.geodesic_points[J],
        ce.source_velocity[J],
        [i.aux[1] for i in ce.geodesic_points[J]],
    )
end

function RadialDiscProfile(ε, ce::CoronaGeodesics; kwargs...)
    J = sortperm(ce.geodesic_points; by = i -> i.x[2])
    radii = @views [i.x[2] for i in ce.geodesic_points[J]]
    times = @views [i.x[1] for i in ce.geodesic_points[J]]

    t = DataInterpolations.LinearInterpolation(times, radii)

    function _emissivity_wrapper(gp)
        ε(gp.x[2])
    end

    function _delay_wrapper(gp)
        t(gp.x[2]) + gp.x[1]
    end

    RadialDiscProfile(_emissivity_wrapper, _delay_wrapper)
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
function emissivity_profile(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    grid = GeometricGrid(),
    N = 100,
    kwargs...,
)
    RadialDiscProfile(tracecorona(m, d, model; kwargs...); grid = grid, N = N)
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
    N = 100,
    λ_max = 10_000.0,
    callback = domain_upper_hemisphere(),
    kwargs...,
)
    δs = deg2rad.(range(0.1, 179.9, N))
    # we assume a point source
    x, v = Gradus.sample_position_velocity(m, model)
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
