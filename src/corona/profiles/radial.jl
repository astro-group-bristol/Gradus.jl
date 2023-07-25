struct RadialDiscProfile{F,R} <: AbstractDiscProfile
    # geodesic point to flux
    f::F
    # geodesic point to time
    t::R
end

function RadialDiscProfile(
    m::AbstractMetric,
    model::AbstractCoronaModel,
    gps::AbstractVector{<:GeodesicPoint},
    source_vels::AbstractVector,
    intensities = nothing;
    kwargs...,
)
    radii = map(i -> i.x[2], gps)
    # ensure sorted: let the user sort so that everything is sure to be
    # in order
    if !issorted(radii)
        error(
            "geodesic points (and therefore source velocities) must be sorted by radii: use `sortperm(gps; by = i -> i.x[2])` to get the sorting permutation for both",
        )
    end

    times = map(i -> i.x[1], gps)
    gs = energy_ratio.(m, gps, source_vels)

    RadialDiscProfile(m, radii, times, gs, intensities; kwargs...)
end

function RadialDiscProfile(rs::AbstractArray, ts::AbstractArray, εs::AbstractArray)
    # create interpolations
    t = DataInterpolations.LinearInterpolation(ts, rs)
    ε = DataInterpolations.LinearInterpolation(εs, rs)
    # wrap geodesic point wrappers
    RadialDiscProfile(gp -> ε(gp.x[2]), gp -> t(gp.x[2]) + gp.x[1])
end

function RadialDiscProfile(rdp::RadialDiscProfile)
    Base.depwarn(
        "This function is deprecated. Note that `emissivity_profile` now returns a `RadialDiscProfile`.",
        :RadialDiscProfile,
    )
    rdp
end

"""
The relativistic correction is calculated via

```math
\\tilde{A} = A \\sqrt{g_{\\mu,\\nu}(x)}
```
"""
function RadialDiscProfile(
    m::AbstractMetric,
    radii::AbstractArray{T},
    times,
    gs,
    intensity,
    ;
    grid = GeometricGrid(),
    N = 100,
) where {T}
    bins = collect(grid(extrema(radii)..., N))

    ibucket = Buckets.IndexBucket(Int, size(radii))
    bucket!(ibucket, Buckets.Simple(), radii, bins)
    groups = Buckets.unpack_bucket(ibucket)

    gs = [@views(mean(gs[g])) for g in groups]
    ts = [@views(mean(times[g])) for g in groups]
    I = if isnothing(intensity)
        # count number of photons in each radial bin
        @. convert(T, length(groups))
    else
        # sum the intensity over the bin
        [@views(sum(intensity[g])) for g in groups]
    end

    # interpolate the energy ratio over the disc
    bin_domain = bins[1:end-1]
    eratios = DataInterpolations.LinearInterpolation(gs, bin_domain)

    for i in eachindex(I)
        R = bins[i]
        r = i == 1 ? 0 : bins[i-1]

        # 2π comes from integrating over dϕ
        dr = R - r
        x = SVector(0, R, π / 2, 0)
        A = dr * _proper_area(m, x)
        # I now stores emissivity
        I[i] = source_to_disc_emissivity(m, I[i], A, x, eratios(R))
    end

    RadialDiscProfile(bins, ts, I)
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

emitted_flux(profile::RadialDiscProfile, gps) = map(profile.f, gps)
delay(profile::RadialDiscProfile, gps) = map(profile.t, gps)
get_emissivity(prof::RadialDiscProfile) = (prof.f.ε.t, prof.f.ε.u)
