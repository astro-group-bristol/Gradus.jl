struct RadialDiscProfile{T,I} <: AbstractDiscProfile
    radii::Vector{T}
    ε::Vector{T}
    t::Vector{T}
    interp_ε::I
    interp_t::I
end

@inline function _enforce_interpolation_bounds(r::Number, prof::RadialDiscProfile)
    r_min = first(prof.radii)
    r_max = last(prof.radii)
    return _enforce_interpolation_bounds(r, r_min, r_max)
end

function emissivity_at(prof::RadialDiscProfile, r::Number)
    r_bounded = _enforce_interpolation_bounds(r, prof)
    prof.interp_ε(r_bounded)
end
emissivity_at(prof::RadialDiscProfile, gp::AbstractGeodesicPoint) =
    emissivity_at(prof, _equatorial_project(gp.x))

function coordtime_at(prof::RadialDiscProfile, r::Number)
    r_bounded = _enforce_interpolation_bounds(r, prof)
    prof.interp_t(r_bounded)
end
coordtime_at(prof::RadialDiscProfile, gp::AbstractGeodesicPoint) =
    coordtime_at(prof, _equatorial_project(gp.x)) + gp.x[1]

function RadialDiscProfile(r, t, ε)
    interp_t = _make_interpolation(view(r, :), view(t, :))
    interp_ε = _make_interpolation(view(r, :), view(ε, :))
    RadialDiscProfile(r, ε, t, interp_ε, interp_t)
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

function RadialDiscProfile(cg::CoronaGeodesics, spec::AbstractCoronalSpectrum; kwargs...)
    J = sortperm(cg.geodesic_points; by = i -> _equatorial_project(i.x))
    @views RadialDiscProfile(
        cg.metric,
        cg.model,
        spec,
        cg.geodesic_points[J],
        cg.source_velocity[J];
        kwargs...,
    )
end

function RadialDiscProfile(
    cg::CoronaGeodesics{<:TraceRadiativeTransfer},
    spec::AbstractCoronalSpectrum;
    kwargs...,
)
    J = sortperm(cg.geodesic_points; by = i -> _equatorial_project(i.x))
    @views RadialDiscProfile(
        cg.metric,
        cg.model,
        spec,
        cg.geodesic_points[J],
        cg.source_velocity[J];
        intensity = [i.aux[1] for i in cg.geodesic_points[J]],
        kwargs...,
    )
end

"""
    TimeDependentRadialDiscProfile{T,I}

Time dependent radial disc profile. Each entry in the radii, emissivity, and
time function maps to a specific contribution, along with the appropriate
weighting when summing the emissivities at similar times together.
"""
struct TimeDependentRadialDiscProfile{T,I} <: AbstractDiscProfile
    weights::Vector{T}
    radii::Vector{Vector{T}}
    ε::Vector{Vector{T}}
    t::Vector{Vector{T}}
end

function emissivity_at(prof::TimeDependentRadialDiscProfile, ρ, t)
    # interpolate the time curve at ρ until we have some (ρ, t)
end
