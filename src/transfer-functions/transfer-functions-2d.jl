function bin_transfer_function(
    time_delays,
    energy,
    flux;
    N_E = 300,
    N_t = 300,
    energy_lims = extrema(energy),
    time_lims = extrema(time_delays),
)
    energy_bins = range(energy_lims..., N_E)
    time_bins = range(time_lims..., N_t)

    de = step(energy_bins)
    dt = step(time_bins)

    transfer_function =
        bucket(energy, time_delays, flux, energy_bins, time_bins; reduction = sum)

    @. transfer_function = transfer_function / (de * dt)
    transfer_function[transfer_function.==0.0] .= NaN
    time_bins, energy_bins, transfer_function
end

function bin_and_interpolate(
    X,
    y::AbstractArray{T};
    log_bins = false,
    nbins = 1000,
    reduction = mean,
) where {T}
    bins = if log_bins
        10 .^ range(log10(minimum(X)), log10(maximum(X)), nbins)
    else
        range(minimum(X), maximum(X), nbins)
    end

    y_binned = bucket(X, y, bins; reduction = reduction)
    DataInterpolations.LinearInterpolation(y_binned, bins)
end

function source_to_disc(
    m::AbstractMetricParams,
    model::AbstractCoronaModel,
    d::AbstractAccretionGeometry,
    max_t;
    n_samples = 10_000,
    sampler = EvenSampler(domain = BothHemispheres(), generator = RandomGenerator()),
    solver_opts...,
)
    sols = tracegeodesics(
        m,
        model,
        d,
        (0.0, max_t);
        n_samples = n_samples,
        save_on = false,
        sampler = sampler,
        solver_opts...,
    )
    points = filter(
        i -> i.status == StatusCodes.IntersectedWithGeometry,
        process_solution.(m, sols.u),
    )
    # sort by radius
    sort!(points, by = i -> i.u2[2])
    points
end

function observer_to_disc(
    m::AbstractMetricParams,
    u,
    plane::AbstractImagePlane,
    d::AbstractAccretionGeometry,
    max_t,
    ;
    solver_opts...,
)
    sols = tracegeodesics(m, u, plane, d, (0.0, max_t); save_on = false, solver_opts...)
    points = process_solution.(m, sols.u)
    points
end

struct LagTransferFunction{T,M,U,P}
    max_t::T
    metric::M
    u::U
    area_elements::Vector{T}
    source_to_disc::Vector{P}
    observer_to_disc::Vector{P}
end

function lagtransfer(
    model::AbstractCoronaModel,
    m::AbstractMetricParams,
    u,
    plane::AbstractImagePlane,
    d::AbstractAccretionGeometry;
    max_t = 2 * u[2],
    n_samples = 10_000,
    sampler = EvenSampler(domain = BothHemispheres(), generator = RandomGenerator()),
    solver_opts...,
)
    s_to_d = source_to_disc(
        m,
        model,
        d,
        max_t;
        n_samples = n_samples,
        sampler = sampler,
        solver_opts...,
    )
    o_to_d = observer_to_disc(
        m,
        u,
        plane,
        d,
        max_t;
        effective_infinity = 1.1 * u[2],
        solver_opts...,
    )

    I = [i.status == StatusCodes.IntersectedWithGeometry for i in o_to_d]
    areas = unnormalized_areas(plane)[I]

    LagTransferFunction(max_t, m, u, areas, s_to_d, o_to_d[I])
end

function _interpolate_profile(tf::LagTransferFunction)
    times = map(i -> i.u2[1], tf.source_to_disc)
    radii = map(i -> i.u2[2], tf.source_to_disc)
    I = sortperm(radii)
    @show any(isnan.(times))
    @show any(isnan.(radii))
    DataInterpolations.LinearInterpolation(times, radii)
end

function binflux(
    tf::LagTransferFunction;
    redshift = ConstPointFunctions.redshift(tf.metric, tf.u),
    kwargs...,
)
    intp = _interpolate_profile(tf)
    r = map(i -> i.u2[2], tf.observer_to_disc)
    t1 = intp.(r)
    t2 = map(i -> i.u2[1], tf.observer_to_disc)
    # add times to get total time
    t = @. t1 + t2

    g = redshift.(tf.metric, tf.observer_to_disc, tf.max_t)

    # calculate flux
    f = @. g^4 * r^(-3) * tf.area_elements / (g * 6.4)
    # normalize
    F = f ./ sum(f)

    tb, eb, td = bin_transfer_function(t, g * 6.4, F; kwargs...)
    # correct for the off by one error in the binning
    tb = tb[1:end-1]
    eb = eb[1:end-1]
    td = td[2:end, 2:end]

    tb, eb, td
end

export bin_transfer_function, lagtransfer, binflux
