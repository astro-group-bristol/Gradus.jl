function render_configuration(
    m,
    position,
    args...;
    image_width,
    image_height,
    fov,
    μ = 0.0,
    q = 0.0,
    trace = TraceGeodesic(; μ = μ, q = q),
    kwargs...,
)
    config, solver_opts = tracing_configuration(
        trace,
        m,
        position,
        _render_velocity_function(m, position, image_width, image_height, fov),
        args...;
        trajectories = image_width * image_height,
        # default do not save full path
        save_on = false,
        kwargs...,
    )
    trace, config, solver_opts
end

function rendergeodesics(
    m::AbstractMetric{T},
    position,
    args...;
    image_width = 350,
    image_height = 250,
    fov = 30.0,
    ensemble = EnsembleEndpointThreads(),
    kwargs...,
) where {T}
    trace, config, solver_opts = render_configuration(
        m,
        position,
        args...;
        ensemble = ensemble,
        image_width,
        image_height,
        fov,
        kwargs...,
    )
    image = zeros(T, (image_height, image_width))
    render_into_image!(image, trace, config; solver_opts...)
    α, β = impact_axes(image_width, image_height, fov)
    α, β, image
end

function prerendergeodesics(
    m::AbstractMetric,
    args...;
    cache = EndpointCache(),
    image_width = 350,
    image_height = 250,
    fov = 3.0,
    kwargs...,
)
    trace, config, solver_opts = render_configuration(
        m,
        position,
        args...;
        image_width,
        image_height,
        fov,
        kwargs...,
    )
    __prerendergeodesics(
        trace,
        config,
        cache;
        image_height = image_height,
        image_width = image_width,
        solver_opts...,
    )
end

function render_into_image!(
    image,
    trace::AbstractTraceParameters,
    config::TracingConfiguration;
    pf = PointFunction((m, gp, λ_max) -> gp.λ_max) ∘
         FilterPointFunction((m, gp, λ_max; kwargs...) -> gp.λ_max < λ_max, NaN),
    solver_opts...,
)
    sol_or_points = __render_geodesics(trace, config; solver_opts...)
    points = sol_or_points_to_points(sol_or_points)
    apply_to_image!(config.metric, image, points, pf, config.λ_domain[2])
    image
end

function apply_to_image!(m::AbstractMetric, image, points, pf, max_time)
    @inbounds Threads.@threads for i in eachindex(points)
        image[i] = pf(m, points[i], max_time)
    end
end

function __prerendergeodesics(
    trace::AbstractTraceParameters,
    config::TracingConfiguration,
    ::SolutionCache;
    image_height,
    image_width,
    kwargs...,
)
    simsols = __render_geodesics(trace, config; kwargs...)
    SolutionRenderCache(config.m, config.λ_domain[2], image_height, image_width, simsols.u)
end

function __prerendergeodesics(
    trace::AbstractTraceParameters,
    config::TracingConfiguration,
    ::EndpointCache;
    image_height,
    image_width,
    kwargs...,
)
    sol_or_points = __render_geodesics(trace, config; kwargs...)
    points = sol_or_points_to_points(sol_or_points)
    EndpointRenderCache(
        config.m,
        config.λ_domain[2],
        image_height,
        image_width,
        reshape(points, (image_height, image_width)),
    )
end

function _render_velocity_function(
    m::AbstractMetric,
    position,
    image_width,
    image_height,
    fov,
)
    y_mid = image_height ÷ 2
    x_mid = image_width ÷ 2
    function velfunc(i)
        Y = i % image_height
        X = i ÷ image_height
        α = x_to_α(X, x_mid, fov)
        β = y_to_β(Y, y_mid, fov)
        map_impact_parameters(m, position, α, β)
    end
end

function __render_geodesics(
    trace::AbstractTraceParameters,
    config::TracingConfiguration;
    verbose = false,
    solver_opts...,
)
    if verbose
        println("+ Starting trace...")
    end

    problem = assemble_tracing_problem(trace, config)
    progress_bar = init_progress_bar("Rendering:", config.trajectories, verbose)
    sol_or_points =
        solve_tracing_problem(problem, config; progress_bar = progress_bar, solver_opts...)

    if verbose
        println("+ Trace complete.")
    end

    sol_or_points
end

sol_or_points_to_points(points::AbstractArray{<:AbstractGeodesicPoint}) = points

sol_or_points_to_points(sols::AbstractArray) = map(process_solution, sols)

sol_or_points_to_points(sol::SciMLBase.EnsembleSolution) = sol_or_points_to_points(sol.u)

export rendergeodesics, prerendergeodesics
