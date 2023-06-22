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

function observer_to_disc(
    m::AbstractMetric,
    u,
    plane::AbstractImagePlane,
    d::AbstractAccretionGeometry,
    max_t,
    ;
    solver_opts...,
)
    points = tracegeodesics(
        m,
        u,
        plane,
        d,
        (0.0, max_t);
        save_on = false,
        ensemble = EnsembleEndpointThreads(),
        solver_opts...,
    )
    points
end

function lagtransfer(
    m::AbstractMetric,
    u,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    plane = PolarPlane(GeometricGrid(); Nr = 800, Nθ = 800, r_max = 50.0),
    max_t = 2 * u[2],
    n_samples = 10_000,
    sampler = EvenSampler(domain = BothHemispheres(), generator = RandomGenerator()),
    verbose = false,
    solver_opts...,
)
    progress_bar = init_progress_bar(
        "Source to disc:   ",
        n_samples + trajectory_count(plane),
        verbose,
    )

    ce = emissivity_profile(
        m,
        d,
        model;
        λ_max = max_t,
        n_samples = n_samples,
        verbose = verbose,
        progress_bar = progress_bar,
        callback = domain_upper_hemisphere(),
        sampler = sampler,
        solver_opts...,
    )

    if verbose
        progress_bar.desc = "Observer to disc: "
    end

    o_to_d = observer_to_disc(
        m,
        u,
        plane,
        d,
        max_t;
        chart = chart_for_metric(m, 1.1 * u[2]),
        verbose = verbose,
        progress_bar = progress_bar,
        callback = domain_upper_hemisphere(),
        solver_opts...,
    )

    # get the indices of all geodesics which intersected with the disc
    I = [i.status == StatusCodes.IntersectedWithGeometry for i in o_to_d]
    areas = unnormalized_areas(plane)[I]

    LagTransferFunction(max_t, u, areas, ce, o_to_d[I])
end

function binflux(tf::LagTransferFunction; kwargs...)
    profile = RadialDiscProfile(r -> r^-3, tf.emissivity_profile)
    binflux(tf, profile; kwargs...)
end


function binflux(
    tf::LagTransferFunction,
    profile::AbstractDiscProfile;
    redshift = ConstPointFunctions.redshift(tf.emissivity_profile.metric, tf.x),
    E₀ = 6.4,
    kwargs...,
)
    t, i_em = delay_flux(profile, tf.observer_to_disc)
    @show t
    g = redshift.(tf.emissivity_profile.metric, tf.observer_to_disc, tf.max_t)

    # calculate flux
    f = @. g^3 * i_em * tf.image_plane_areas
    # normalize
    F = f ./ sum(f)

    tb, eb, td = bin_transfer_function(t, g * E₀, F; kwargs...)
    # subtract initial time
    tb .- tf.x[2], eb, td
end

export bin_transfer_function, lagtransfer, binflux
