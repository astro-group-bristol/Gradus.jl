"""
    source_to_disc_emissivity(m, 𝒩, A, x, g)

Compute the emissivity (in arbitrary units) in the disc element with area `A`, photon
count `𝒩`, central position `x`, and redshift `g`. Evaluates

```math
\\varepsilon = \\frac{\\mathscr{N}}{\\tilde{A} g^2},
```

where ``\\tilde{A}`` is the relativistically corrected area of `A`. The relativistic correction
is calculated via

```math
\\tilde{A} = A \\sqrt{g_{\\mu,\\nu}(x)}
```
"""
function source_to_disc_emissivity(m::AbstractStaticAxisSymmetric, 𝒩, A, x, g)
    gcomp = metric_components(m, x)
    # account for relativistic effects in area
    A_corrected = A * √(gcomp[2] * gcomp[3])
    # divide by area to get number density
    𝒩 / (g^2 * A_corrected)
end

struct CoronalEmissivity{T,M,G,C,P,V}
    trace::T
    metric::M
    geometry::G
    model::C
    geodesic_points::P
    source_velocity::V
end

function emissivity_profile(
    m::AbstractMetric,
    g::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    λ_max = 10_000,
    n_samples = 1024,
    sampler = EvenSampler(domain = BothHemispheres(), generator = RandomGenerator()),
    trace = TraceGeodesic(),
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
        callback = domain_upper_hemisphere(),
        kwargs...,
    )
    mask = [i.status == StatusCodes.IntersectedWithGeometry for i in gps]
    CoronalEmissivity(trace, m, g, model, gps[mask], source_vels[mask])
end

function RadialDiscProfile(ce::CoronalEmissivity; kwargs...)
    J = sortperm(ce.geodesic_points; by = i -> i.x[2])
    @views RadialDiscProfile(
        ce.metric,
        ce.model,
        ce.geodesic_points[J],
        ce.source_velocity[J],
    )
end

function RadialDiscProfile(ce::CoronalEmissivity{<:TraceRadiativeTransfer}; kwargs...)
    J = sortperm(ce.geodesic_points; by = i -> i.x[2])
    @views RadialDiscProfile(
        ce.metric,
        ce.model,
        ce.geodesic_points[J],
        ce.source_velocity[J],
        [i.aux[1] for i in ce.geodesic_points[J]],
    )
end

function RadialDiscProfile(ε, ce::CoronalEmissivity; kwargs...)
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
