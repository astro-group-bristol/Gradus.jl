function __rendergeodesics(
    m::AbstractMetricParams{T},
    init_pos;
    image_width,
    image_height,
    kwargs...,
) where {T}
    image = zeros(T, (image_height, image_width))
    render_into_image!(
        m,
        init_pos,
        image;
        image_width = image_width,
        image_height = image_height,
        kwargs...,
    )
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
    cache::DefaultCache;
    kwargs...,
) where {T}
    # default is endpoints for now
    __prerendergeodesics(m, init_pos, EndpointCache(); kwargs...)
end

function __prerendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    cache::SolutionCache;
    kwargs...,
) where {T}
    prerender_solutions(m, init_pos; kwargs...)
end

function __prerendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    cache::EndpointCache;
    kwargs...,
) where {T}
    prerender_endpoints(m, init_pos; kwargs...)
end

#get_parameter_type(m, u, v, max_time) =
#    typeof(GeodesicTracer.integrator_problem(m, u, v, (0.0, max_time)).p)

function prerender_endpoints(
    m::AbstractMetricParams{T},
    init_pos::V;
    image_width,
    image_height,
    fov_factor,
    max_time,
    solver_opts...,
) where {T,V}
    y_mid = image_height ÷ 2
    x_mid = image_width ÷ 2

    #u_cache = Matrix{V}(undef, (image_height, image_width))
    #v_cache = Matrix{V}(undef, (image_height, image_width))
    #t_cache = Matrix{T}(undef, (image_height, image_width))
    #p_cache = Matrix{get_parameter_type(m, init_pos, init_pos, max_time)}(undef, (image_height, image_width))
    #retcode_cache = Matrix{Symbol}(undef, (image_height, image_width))

    point_cache =
        Matrix{GradusBase.geodesic_point_type(m)}(undef, (image_height, image_width))

    vs = fill(init_pos, image_width)
    us = fill(init_pos, image_width)

    for Y = 1:image_height

        β = T(y_to_β(Y, y_mid, fov_factor))
        α_generator_row = (T(x_to_α(X, x_mid, fov_factor)) for X = 1:image_width)

        calculate_velocities!(vs, m, init_pos, α_generator_row, β)
        simsols =
            tracegeodesics(m, us, vs, (T(0.0), max_time); save_on = false, solver_opts...)
        println("+ $Y / $image_height ...")

        point_cache[Y, :] = get_endpoint(m, simsols)
    end

    EndpointRenderCache(m, max_time, image_height, image_width, point_cache)
end

function prerender_solutions(
    m::AbstractMetricParams{T},
    init_pos;
    image_width,
    image_height,
    fov_factor,
    max_time,
    solver_opts...,
) where {T}
    y_mid = image_height ÷ 2
    x_mid = image_width ÷ 2


    # this hits the garbage collector like a truck
    # need a better pre-cache for the velocity and position vectors
    # maybe using a row by row renderer as before? 
    # or can we somehow repeat us without copying, as with FillArrays.jl?
    # -- tried that, doesn't make enough of a difference

    vs = fill(init_pos, image_width)
    us = fill(init_pos, image_width)

    simsol_array = map(1:image_height) do Y

        β = T(y_to_β(Y, y_mid, fov_factor))
        α_generator_row = (T(x_to_α(X, x_mid, fov_factor)) for X = 1:image_width)

        calculate_velocities!(vs, m, init_pos, α_generator_row, β)
        simsols = tracegeodesics(
            m,
            us,
            vs,
            (T(0.0), max_time);
            save_on = false,
            solver_opts...,
        )
        println("+ $Y / $image_height ...")
        simsols
    end

    SolutionRenderCache(m, max_time, image_height, image_width, simsol_array)
end

function render_into_image!(
    m::AbstractMetricParams{T},
    init_pos,
    image;
    image_width,
    image_height,
    fov_factor,
    max_time,
    pf,
    verbose = true,
    solver_opts...,
) where {T}
    y_mid = image_height ÷ 2
    x_mid = image_width ÷ 2

    vs = fill(init_pos, image_width)
    us = fill(init_pos, image_width)
    for Y = 1:image_height

        β = T(y_to_β(Y, y_mid, fov_factor))
        α_generator_row = (T(x_to_α(X, x_mid, fov_factor)) for X = 1:image_width)

        calculate_velocities!(vs, m, init_pos, α_generator_row, β)

        # do the render
        simsols = tracegeodesics(
            m,
            us,
            vs,
            (T(0.0), max_time);
            save_on = false,
            verbose = verbose,
            solver_opts...,
        )
        apply_to_location!(m, @view(image[Y, :]), simsols, pf, max_time)

        if verbose
            println("+ $Y / $image_height ...")
        end
    end

    image
end

function apply_to_location!(m::AbstractMetricParams{T}, loc, sols, pf, max_time) where {T}
    @threads for i = 1:length(sols)
        loc[i] = pf(m, get_endpoint(m, sols[i]), max_time)
    end
end
