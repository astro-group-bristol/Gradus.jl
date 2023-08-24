struct RadialDiscProfile{F,R} <: AbstractDiscProfile
    # geodesic point to flux
    f::F
    # geodesic point to time
    t::R
end

function RadialDiscProfile(rs::AbstractArray, ts::AbstractArray, εs::AbstractArray)
    # create interpolations
    t = _make_interpolation(rs, ts)
    ε = _make_interpolation(rs, εs)
    # wrap geodesic point wrappers
    RadialDiscProfile(
        gp -> ε(_equatorial_project(gp.x)),
        gp -> t(_equatorial_project(gp.x)) + gp.x[1],
    )
end

_get_grouped_intensity(T::Type, groupings, ::Nothing) = @. convert(T, length(groupings))
_get_grouped_intensity(::Type, groupings, intensity) =
    [@views(sum(intensity[grouping])) for grouping in groupings]

function _build_radial_profile(
    m::AbstractMetric,
    spec::AbstractCoronalSpectrum,
    radii,
    times,
    source_velocities,
    points::AbstractVector{<:AbstractGeodesicPoint{T}},
    intensity;
    grid = GeometricGrid(),
    N = 100,
) where {T}
    @assert size(radii) == size(source_velocities)
    @assert size(points) == size(source_velocities)

    # function for obtaining keplerian velocities
    _disc_velocity = _keplerian_velocity_projector(m)

    bins = grid(extrema(radii)..., N) |> collect

    # find the grouping with an index bucket
    # that is, the points with the same r in Δr on the disc
    ibucket = Buckets.IndexBucket(Int, size(radii), length(bins))
    bucket!(ibucket, Buckets.Simple(), radii, bins)
    groupings = Buckets.unpack_bucket(ibucket)
    @show groupings

    # need to interpolate the redshifts, so calculate those first
    gs = map(groupings) do grouping
        g_total = zero(T)
        for i in grouping
            gp = points[i]
            v_disc = _disc_velocity(gp.x)
            g_total += energy_ratio(m, gp, source_velocities[i], v_disc)
        end
        # get the mean
        g_total / length(grouping)
    end

    g_interp = _make_interpolation(bins, gs)
    grouped_I = _get_grouped_intensity(T, groupings, intensity)

    for (i, I) in enumerate(grouped_I)
        R = bins[i]
        r = i == 1 ? 0 : bins[i-1]
        dr = R - r
        x = SVector(0, R, π / 2, 0)
        v_disc = _disc_velocity(x)
        A = dr * _proper_area(m, x)
        # now stores emissivity
        grouped_I[i] = source_to_disc_emissivity(m, spec, I, A, x, g_interp(R), v_disc)
    end

    ts = [@views(mean(times[grouping])) for grouping in groupings]

    bins, ts, grouped_I
end

function RadialDiscProfile(
    m::AbstractMetric,
    model::AbstractCoronaModel,
    spec::AbstractCoronalSpectrum,
    points::AbstractVector{<:AbstractGeodesicPoint},
    source_velocities::AbstractVector;
    intensity = nothing,
    kwargs...,
)
    radii = map(i -> _equatorial_project(i.x), points)
    # ensure sorted: let the user sort so that everything is sure to be
    # in order
    if !issorted(radii)
        error(
            "geodesic points (and therefore also source velocities) must be sorted by radii: use `sortperm(points; by = i -> i.x[2])` to get the sorting permutation for both",
        )
    end

    times = map(i -> i.x[1], points)
    r, t, ε = _build_radial_profile(
        m,
        spec,
        radii,
        times,
        source_velocities,
        points,
        intensity;
        kwargs...,
    )
    RadialDiscProfile(r, t, ε)
end

function RadialDiscProfile(rdp::RadialDiscProfile)
    Base.depwarn(
        "This function is deprecated. Note that `emissivity_profile` now returns a `RadialDiscProfile`.",
        :RadialDiscProfile,
    )
    rdp
end

function RadialDiscProfile(ce::CoronaGeodesics, spec::AbstractCoronalSpectrum; kwargs...)
    J = sortperm(ce.geodesic_points; by = i -> _equatorial_project(i.x))
    @views RadialDiscProfile(
        ce.metric,
        ce.model,
        spec,
        ce.geodesic_points[J],
        ce.source_velocity[J];
        kwargs...,
    )
end

function RadialDiscProfile(
    ce::CoronaGeodesics{<:TraceRadiativeTransfer},
    spec::AbstractCoronalSpectrum;
    kwargs...,
)
    J = sortperm(ce.geodesic_points; by = i -> _equatorial_project(i.x))
    @views RadialDiscProfile(
        ce.metric,
        ce.model,
        spec,
        ce.geodesic_points[J],
        ce.source_velocity[J];
        intensity = [i.aux[1] for i in ce.geodesic_points[J]],
        kwargs...,
    )
end

emitted_flux(profile::RadialDiscProfile, gps) = map(profile.f, gps)
delay(profile::RadialDiscProfile, gps) = map(profile.t, gps)
get_emissivity(prof::RadialDiscProfile) = (prof.f.ε.t, prof.f.ε.u)
