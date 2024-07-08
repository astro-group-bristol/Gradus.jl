struct RadialDiscProfile{V<:AbstractVector,I} <: AbstractDiscProfile
    radii::V
    ε::V
    t::V
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
struct TimeDependentRadialDiscProfile{T} <: AbstractDiscProfile
    weights::Vector{T}
    radii::Vector{Vector{T}}
    t::Vector{Vector{T}}
    ε::Vector{Vector{T}}
end

# effectively time averaged values
function emissivity_at(prof::TimeDependentRadialDiscProfile{T}, ρ) where {T}
    sum(eachindex(prof.radii)) do i
        radii = prof.radii[i]
        if (ρ >= radii[1]) && (ρ <= radii[end])
            _make_interpolation(radii, prof.ε[i])(ρ)
        else
            zero(T)
        end
    end
end

function emissivity_interp(prof::TimeDependentRadialDiscProfile{T}, ρ) where {T}
    ts = zeros(T, length(prof.weights))
    εs = zeros(T, length(prof.weights))
    # interpolate the time curve at ρ until we have some (ρ, t)
    for i in eachindex(prof.radii)
        radii = prof.radii[i]
        if ρ >= radii[1] && ρ <= radii[end]
            t_interp = _make_interpolation(radii, prof.t[i])
            ε_interp = _make_interpolation(radii, prof.ε[i])
            ts[i] = t_interp(ρ)
            εs[i] = ε_interp(ρ)
        else
            ts[i] = NaN
            εs[i] = NaN
        end
    end

    J = sortperm(ts)
    _make_interpolation(ts[J], εs[J])
end

function emissivity_interp_limits(prof::TimeDependentRadialDiscProfile{T}, ρ) where {T}
    ts = map(eachindex(prof.radii)) do i
        radii = prof.radii[i]
        if ρ >= radii[1] && ρ <= radii[end]
            t_interp = _make_interpolation(radii, prof.t[i])
            t_interp(ρ)
        else
            NaN
        end
    end
    filter!(!isnan, ts)
    if length(ts) > 0
        extrema(ts)
    else
        (zero(T), zero(T))
    end
end

"""
    struct RingCoronaProfile{T} <: AbstractDiscProfile

A specialised disc profile that can stores the various
[`TimeDependentRadialDiscProfile`](@ref) for the ring-like extended corona.
"""
struct RingCoronaProfile{T} <: AbstractDiscProfile
    left_arm::TimeDependentRadialDiscProfile{T}
    right_arm::TimeDependentRadialDiscProfile{T}
end

function emissivity_at(prof::RingCoronaProfile{T}, ρ) where {T}
    emissivity_at(prof.left_arm, ρ) + emissivity_at(prof.right_arm, ρ)
end

function emissivity_interp(prof::RingCoronaProfile{T}, ρ) where {T}
    left_arm = emissivity_interp(prof.left_arm, ρ)
    right_arm = emissivity_interp(prof.right_arm, ρ)
    function _add_arms(t)
        left = if t >= left_arm.t[1] && t <= left_arm.t[end]
            left_arm(t)
        else
            zero(T)
        end
        right = if t >= right_arm.t[1] && t <= right_arm.t[end]
            right_arm(t)
        else
            zero(T)
        end

        left + right
    end
end

function emissivity_interp_limits(prof::RingCoronaProfile, ρ)
    l_min, l_max = emissivity_interp_limits(prof.left_arm, ρ)
    r_min, r_max = emissivity_interp_limits(prof.right_arm, ρ)
    min(l_min, r_min), max(l_max, r_max)
end

struct DiscCoronaProfile{T} <: AbstractDiscProfile
    radii::Vector{T}
    rings::Vector{RingCoronaProfile{T}}
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(prof::DiscCoronaProfile))
    print(
        io,
        """DiscCoronaProfile
  . N rings      : $(length(prof.radii))
  . r (min, max) : $(extrema(prof.radii))
""",
    )
end

function _ring_weighting(prof::DiscCoronaProfile, i)
    δr = prof.radii[2] - prof.radii[1]
    prof.radii[i] * δr
end

function emissivity_at(prof::DiscCoronaProfile{T}, ρ) where {T}
    sum(eachindex(prof.radii)) do i
        emissivity_at(prof.rings[i], ρ) * _ring_weighting(prof, i)
    end
end

function emissivity_interp(prof::DiscCoronaProfile{T}, ρ) where {T}
    funcs = [emissivity_interp(ring, ρ) for ring in prof.rings]
    function _sum_ring_contributions(t)
        sum(eachindex(funcs)) do i
            funcs[i](t) * _ring_weighting(prof, i)
        end
    end
end

function emissivity_interp_limits(prof::DiscCoronaProfile, ρ)
    _min, _max = emissivity_interp_limits(prof.rings[1], ρ)
    for i = 2:lastindex(prof.rings)
        l_min, l_max = emissivity_interp_limits(prof.rings[i], ρ)
        _min = min(_min, l_min)
        _max = max(_max, l_max)
    end
    (_min, _max)
end
