function rendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    max_time;
    pf = PointFunction((m, gp, mt) -> gp.t2) ∘
         FilterPointFunction((m, gp, max_time; kwargs...) -> gp.t2 < max_time, NaN),
    image_width = 350,
    image_height = 250,
    fov_factor = 3.0,
    kwargs...,
) where {T}
    image = zeros(T, (image_height, image_width))
    render_into_image!(
        image,
        m,
        init_pos,
        max_time;
        pf = pf,
        image_width = image_width,
        image_height = image_height,
        fov_factor = fov_factor,
        kwargs...,
    )
    image
end

function prerendergeodesics(
    m::AbstractMetricParams,
    init_pos,
    max_time;
    cache = EndpointCache(),
    image_width = 350,
    image_height = 250,
    fov_factor = 3.0,
    kwargs...,
)
    __prerendergeodesics(
        m,
        init_pos,
        max_time,
        cache;
        image_width = image_width,
        image_height = image_height,
        fov_factor = fov_factor,
        kwargs...,
    )
end

function render_into_image!(
    image,
    m::AbstractMetricParams,
    init_pos,
    max_time;
    pf,
    kwargs...,
)
    sol_or_points = __render_geodesics(m, init_pos, max_time, ; kwargs...)
    points = sol_or_points_to_points(sol_or_points)
    apply_to_image!(m, image, points, pf, max_time)
    image
end

function apply_to_image!(m::AbstractMetricParams, image, points, pf, max_time)
    @inbounds Threads.@threads for i in eachindex(points)
        image[i] = pf(m, points[i], max_time)
    end
end

function __prerendergeodesics(
    m::AbstractMetricParams,
    init_pos,
    cache::AbstractRenderCache;
    kwargs...,
)
    error("Not implemented for render cache strategy '$typeof(cache)'.")
end

function __prerendergeodesics(
    m::AbstractMetricParams,
    init_pos,
    max_time,
    cache::SolutionCache;
    image_height,
    image_width,
    kwargs...,
)
    simsols = __render_geodesics(
        m,
        init_pos,
        max_time;
        image_height = image_height,
        image_width = image_width,
        ensemble = EnsembleThreads(),
        kwargs...,
    )
    SolutionRenderCache(m, max_time, image_height, image_width, simsols.u)
end

function __prerendergeodesics(
    m::AbstractMetricParams,
    init_pos,
    max_time,
    cache::EndpointCache;
    image_height,
    image_width,
    kwargs...,
)
    sol_or_points = __render_geodesics(
        m,
        init_pos,
        max_time;
        image_height = image_height,
        image_width = image_width,
        kwargs...,
    )

    points = sol_or_points_to_points(sol_or_points)

    EndpointRenderCache(
        m,
        max_time,
        image_height,
        image_width,
        reshape(points, (image_height, image_width)),
    )
end

function __render_geodesics(
    m::AbstractMetricParams{T},
    init_pos::AbstractVector{T},
    max_time;
    image_width,
    image_height,
    fov_factor,
    verbose = false,
    ensemble = EnsembleEndpointThreads(),
    solver_opts...,
) where {T}
    y_mid = image_height ÷ 2
    x_mid = image_width ÷ 2

    trajectories = image_width * image_height

    if verbose
        println("+ Starting trace...")
    end

    progress_bar = init_progress_bar("Rendering:", trajectories, verbose)

    function velfunc(i)
        Y = i % image_height
        X = i ÷ image_height
        α = x_to_α(X, x_mid, fov_factor)
        β = y_to_β(Y, y_mid, fov_factor)
        map_impact_parameters(m, init_pos, α, β)
    end

    sol_or_points = tracegeodesics(
        m,
        init_pos,
        velfunc,
        (T(0.0), max_time);
        save_on = false,
        verbose = verbose,
        trajectories = trajectories,
        progress_bar = progress_bar,
        ensemble = ensemble,
        solver_opts...,
    )

    if verbose
        println("+ Trace complete.")
    end

    sol_or_points
end

function sol_or_points_to_points(points::AbstractArray{<:AbstractGeodesicPoint})
    points
end

function sol_or_points_to_points(sols::AbstractArray)
    map(process_solution, sols)
end

function sol_or_points_to_points(sol::SciMLBase.EnsembleSolution)
    sol_or_points_to_points(sol.u)
end

export rendergeodesics, prerendergeodesics
