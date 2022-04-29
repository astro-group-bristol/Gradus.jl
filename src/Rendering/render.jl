function render_into_image!(
    image,
    m::AbstractMetricParams{T},
    init_pos,
    max_time;
    pf,
    kwargs...,
) where {T}
    simsols = __render_geodesics(m, init_pos, max_time, ; kwargs...)
    apply_to_image!(m, image, simsols, pf, max_time)
    image
end

function apply_to_image!(m::AbstractMetricParams{T}, image, sols, pf, max_time) where {T}
    @inbounds @threads for i = 1:length(sols)
        image[i] = pf(m, getendpoint(m, sols[i]), max_time)
    end
end

function __prerendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    cache::AbstractRenderCache;
    kwargs...,
) where {T}
    error("Not implemented for render cache strategy '$typeof(cache)'.")
end

function __prerendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    max_time,
    cache::SolutionCache;
    image_height,
    image_width,
    kwargs...,
) where {T}
    simsols = __render_geodesics(
        m,
        init_pos,
        max_time;
        image_height = image_height,
        image_width = image_width,
        kwargs...,
    )
    SolutionRenderCache(m, max_time, image_height, image_width, simsols.u)
end

function __prerendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    max_time,
    cache::EndpointCache;
    image_height,
    image_width,
    kwargs...,
) where {T}
    simsols = __render_geodesics(
        m,
        init_pos,
        max_time;
        image_height = image_height,
        image_width = image_width,
        kwargs...,
    )

    point_cache = getendpoint(m, simsols)

    EndpointRenderCache(
        m,
        max_time,
        image_height,
        image_width,
        reshape(point_cache, (image_height, image_width)),
    )
end

function __render_geodesics(
    m::AbstractMetricParams{T},
    init_pos,
    max_time;
    image_width,
    image_height,
    fov_factor,
    verbose = false,
    solver_opts...,
) where {T}
    y_mid = image_height ÷ 2
    x_mid = image_width ÷ 2

    trajectories = image_width * image_height

    if verbose
        println("+ Starting trace...")
    end

    progress_bar = ProgressMeter.Progress(
        trajectories;
        desc = "Rendering:",
        dt = 0.5,
        barglyphs = BarGlyphs("[=> ]"),
        barlen = 40,
        color = :none,
        enabled = verbose,
    )

    function velfunc(i)
        Y = i % image_height
        X = i ÷ image_height
        α = x_to_α(X, x_mid, fov_factor)
        β = y_to_β(Y, y_mid, fov_factor)
        ProgressMeter.next!(progress_bar)
        map_impact_parameters(m, init_pos, α, β)
    end

    simsols = tracegeodesics(
        m,
        init_pos,
        velfunc,
        (T(0.0), max_time);
        save_on = false,
        verbose = verbose,
        trajectories = trajectories,
        solver_opts...,
    )

    if verbose
        println("+ Trace complete.")
    end

    simsols
end
