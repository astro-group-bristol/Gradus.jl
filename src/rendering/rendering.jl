function render_configuration(
    m::AbstractMetric{T},
    position,
    args...;
    image_width,
    image_height,
    αlims,
    βlims,
    μ = 0,
    q = 0,
    trace = TraceGeodesic(; μ = μ, q = q),
    kwargs...,
) where {T}
    config, solver_opts = tracing_configuration(
        trace,
        m,
        position,
        _render_velocity_function(m, position, image_width, image_height, αlims, βlims),
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
    image_width = 375,
    image_height = 250,
    αlims = (-60, 60),
    βlims = (-40, 40),
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
        αlims,
        βlims,
        kwargs...,
    )
    image = zeros(T, (image_height, image_width))
    render_into_image!(image, trace, config; solver_opts...)
    α, β = impact_axes(image_width, image_height, αlims, βlims)
    α, β, image
end

function prerendergeodesics(
    m::AbstractMetric{T},
    position,
    args...;
    cache = EndpointCache(),
    image_width = 375,
    image_height = 250,
    αlims = (-60, 60),
    βlims = (-40, 40),
    kwargs...,
) where {T}
    trace, config, solver_opts = render_configuration(
        m,
        position,
        args...;
        image_width,
        image_height,
        αlims,
        βlims,
        kwargs...,
    )
    cache = _prerendergeodesics(
        trace,
        config,
        cache;
        image_height = image_height,
        image_width = image_width,
        solver_opts...,
    )
    α, β = impact_axes(image_width, image_height, αlims, βlims)
    α, β, cache
end

function render_into_image!(
    image,
    trace::AbstractTrace,
    config::TracingConfiguration{T};
    pf = PointFunction((m, gp, λ_max) -> gp.λ_max) ∘
         FilterPointFunction((m, gp, λ_max; kwargs...) -> gp.λ_max < λ_max, T(NaN)),
    solver_opts...,
) where {T}
    sol_or_points = _render_geodesics(trace, config; solver_opts...)
    points = sol_or_points_to_points(sol_or_points)
    apply_to_image!(config.metric, image, points, pf, config.λ_domain[2])
    image
end

function apply_to_image!(m::AbstractMetric, image, points, pf, max_time)
    @inbounds Threads.@threads for i in eachindex(points)
        image[i] = pf(m, points[i], max_time)
    end
end

function _prerendergeodesics(
    trace::AbstractTrace,
    config::TracingConfiguration,
    ::SolutionCache;
    image_height,
    image_width,
    kwargs...,
)
    simsols = _render_geodesics(trace, config; kwargs...)
    SolutionRenderCache(config.m, config.λ_domain[2], image_height, image_width, simsols.u)
end

function _prerendergeodesics(
    trace::AbstractTrace,
    config::TracingConfiguration,
    ::EndpointCache;
    image_height,
    image_width,
    kwargs...,
)
    sol_or_points = _render_geodesics(trace, config; kwargs...)
    points = sol_or_points_to_points(sol_or_points)
    EndpointRenderCache(
        config.metric,
        config.λ_domain[2],
        image_height,
        image_width,
        reshape(points, (image_height, image_width)),
    )
end

function _render_velocity_function(
    m::AbstractMetric{T},
    position,
    image_width,
    image_height,
    αlims,
    βlims,
) where {T}
    @assert issorted(αlims) "α limits must be sorted"
    @assert issorted(βlims) "β limits must be sorted"

    αs = range(T(αlims[1]), T(αlims[2]), image_width)
    βs = range(T(βlims[1]), T(βlims[2]), image_height)
    xfm = lnr_momentum_to_global_velocity_transform(m, position)
    function velfunc(i)
        # get index on image plane
        x = (i - 1) ÷ image_height + 1
        y = mod1(i, image_height)
        # offset a little to avoid coordinate singularities when α = 0
        α = αs[x] + T(1e-6)
        β = βs[y] + T(1e-6)
        xfm(local_momentum(position[2], α, β))
    end
end

function _render_geodesics(
    trace::AbstractTrace,
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

sol_or_points_to_points(sols::AbstractArray) = map(unpack_solution, sols)

sol_or_points_to_points(sol::SciMLBase.EnsembleSolution) = sol_or_points_to_points(sol.u)

export rendergeodesics, prerendergeodesics
