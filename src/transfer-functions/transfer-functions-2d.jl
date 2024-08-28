function InterpolatingTransferBranches(branches)
    radii = map(branch -> branch.rₑ, branches)
    if !issorted(radii)
        I = sortperm(radii)
        radii = radii[I]
        _branches = branches[I]
    else
        _branches = branches
    end

    gmin = map(i -> i.gmin, _branches)
    gmax = map(i -> i.gmax, _branches)
    InterpolatingTransferBranches(_branches, radii, gmin, gmax)
end

function Base.show(io::IO, ::MIME"text/plain", itb::InterpolatingTransferBranches)
    text = """InterpolatingTransferBranches
      . branches   :  $(length(itb.branches))
      . rₑ extrema :  $(extrema(itb.radii))"""
    print(io, text)
end

function _lazy_interpolate(f1, f2, θ)
    function _lazy_interpolate_kernel(index, weight)
        y1 = _linear_interpolate(f1.u, index, weight)
        y2 = _linear_interpolate(f2.u, index, weight)
        _linear_interpolate(y1, y2, θ)
    end
    function _lazy_interpolate_kernel(x)
        _linear_interpolate(f1(x), f2(x), θ)
    end
end

function _lazy_interpolate(branch::Vector{<:TransferBranches}, idx, θ)
    b1 = branch[idx]
    b2 = branch[idx+1]
    (
        _lazy_interpolate(b1.lower_f, b2.lower_f, θ),
        _lazy_interpolate(b1.upper_f, b2.upper_f, θ),
        _lazy_interpolate(b1.lower_t, b2.lower_t, θ),
        _lazy_interpolate(b1.upper_t, b2.upper_t, θ),
    )
end

function _lazy_interpolate(grid::CunninghamTransferGrid, idx, θ)
    x = grid.g✶_grid
    @views (
        _lazy_interpolate(
            _make_interpolation(x, grid.lower_f[:, idx]),
            _make_interpolation(x, grid.lower_f[:, idx+1]),
            θ,
        ),
        _lazy_interpolate(
            _make_interpolation(x, grid.upper_f[:, idx]),
            _make_interpolation(x, grid.upper_f[:, idx+1]),
            θ,
        ),
        _lazy_interpolate(
            _make_interpolation(x, grid.lower_time[:, idx]),
            _make_interpolation(x, grid.lower_time[:, idx+1]),
            θ,
        ),
        _lazy_interpolate(
            _make_interpolation(x, grid.upper_time[:, idx]),
            _make_interpolation(x, grid.upper_time[:, idx+1]),
            θ,
        ),
    )
end

# interpolate over radial coordinate
function (itb::InterpolatingTransferBranches)(r)
    idx = max(
        1,
        min(
            DataInterpolations.searchsortedlastcorrelated(itb.radii, r, 0),
            length(itb.radii) - 1,
        ),
    )
    r1, r2 = itb.radii[idx], itb.radii[idx+1]
    # interpolation weight
    θ = (r - r1) / (r2 - r1)

    gmin = _linear_interpolate(itb.gmin, idx, θ)
    gmax = _linear_interpolate(itb.gmax, idx, θ)
    upper_f, lower_f, upper_t, lower_t = _lazy_interpolate(itb.branches, idx, θ)
    TransferBranches{false}(upper_f, lower_f, upper_t, lower_t, gmin, gmax, r)
end

function (grid::CunninghamTransferGrid)(r)
    idx = max(
        1,
        min(
            DataInterpolations.searchsortedlastcorrelated(grid.r_grid, r, 0),
            length(grid.r_grid) - 1,
        ),
    )
    r1, r2 = grid.r_grid[idx], grid.r_grid[idx+1]
    # interpolation weight
    θ = (r - r1) / (r2 - r1)

    gmin = _linear_interpolate(grid.g_min, idx, θ)
    gmax = _linear_interpolate(grid.g_max, idx, θ)
    upper_f, lower_f, upper_t, lower_t = _lazy_interpolate(grid, idx, θ)
    TransferBranches{true}(upper_f, lower_f, upper_t, lower_t, gmin, gmax, r)
end

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

    ce = tracecorona(
        m,
        d,
        model;
        λmax = max_t,
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
    I = [status(i) == StatusCodes.IntersectedWithGeometry for i in o_to_d]
    areas = unnormalized_areas(plane)[I]

    LagTransferFunction(max_t, u, areas, ce, o_to_d[I])
end

function binflux(tf::LagTransferFunction; kwargs...)
    profile = AnalyticRadialDiscProfile(r -> r^-3, tf.coronal_geodesics)
    binflux(tf, profile; kwargs...)
end

function binflux(
    tf::LagTransferFunction,
    profile::AbstractDiscProfile;
    redshift = ConstPointFunctions.redshift(tf.coronal_geodesics.metric, tf.x),
    E₀ = 6.4,
    t0 = tf.x[2],
    kwargs...,
)
    t = coordtime_at(profile, tf.observer_to_disc)
    ε = emissivity_at(profile, tf.observer_to_disc)
    g = redshift.(tf.coronal_geodesics.metric, tf.observer_to_disc, tf.max_t)
    # calculate flux
    f = @. g^3 * ε * tf.image_plane_areas
    # normalize
    F = f ./ sum(f)

    tb, eb, td = bin_transfer_function(t, g * E₀, F; kwargs...)
    # subtract initial time
    tb .- t0, eb, td
end

export bin_transfer_function, lagtransfer, binflux
