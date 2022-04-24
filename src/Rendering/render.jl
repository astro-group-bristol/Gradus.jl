function render_into_image!(
    image,
    m::AbstractMetricParams{T},
    init_pos,
    max_time
    ;
    pf,
    kwargs...,
) where {T}
    simsols = __render_geodesics(
        m,
        init_pos,
        max_time,
        ;
        kwargs...
    )
    apply_to_image!(m, image, simsols, pf, max_time)
    image
end

function apply_to_image!(m::AbstractMetricParams{T}, image, sols, pf, max_time) where {T}
    @inbounds @threads for i = 1:length(sols)
        image[i] = pf(m, get_endpoint(m, sols[i]), max_time)
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
    kwargs...,
) where {T}
    simsols = __render_geodesics(m, init_pos, max_time; kwargs...)
    SolutionRenderCache(m, max_time, image_height, image_width, simsols.u)
end

function __prerendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    max_time,
    cache::EndpointCache;
    kwargs...,
) where {T}
simsols = __render_geodesics(m, init_pos, max_time; kwargs...)

    point_cache = get_endpoint(m, simsols)
    reshape!(point_cache, (image_height, image_width))

    EndpointRenderCache(m, max_time, image_height, image_width, point_cache)
end

function __render_geodesics(
    m::AbstractMetricParams{T},
    init_pos,
    max_time
    ;
    image_width,
    image_height,
    fov_factor,
    verbose = true,
    solver_opts...,
) where {T}
    y_mid = image_height ÷ 2
    x_mid = image_width ÷ 2

    function velfunc(i)
        X = i % image_width
        Y = i ÷ image_width
        α = x_to_α(X, x_mid, fov_factor)
        β = y_to_β(Y, y_mid, fov_factor)
        map_impact_parameters(m, init_pos, α, β)
    end

    simsols = tracegeodesics(
        m,
        init_pos,
        velfunc,
        (T(0.0), max_time);
        save_on = false,
        verbose = verbose,
        trajectories = image_width * image_height,
        solver_opts...,
    )
    
    simsols
end
