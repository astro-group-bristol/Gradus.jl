function geodesic_ode_problem(
    ::TraceGeodesic,
    m::AbstractMetric,
    pos,
    vel,
    time_domain,
    callback,
)
    function f(u::SVector{8,T}, p, λ) where {T}
        @inbounds let x = SVector{4,T}(@view(u[1:4])), v = SVector{4,T}(@view(u[5:8]))
            dv = SVector{4,T}(geodesic_equation(m, x, v))
            vcat(v, dv)
        end
    end

    u_init = vcat(pos, vel)
    ODEProblem{false}(
        f,
        u_init,
        time_domain,
        IntegrationParameters(StatusCodes.NoStatus);
        callback = callback,
    )
end

function assemble_tracing_problem(
    trace::AbstractTrace,
    config::TracingConfiguration,
)
    # create the callback set for the problem
    cbs = create_callback_set(config.metric, config.callback, config.chart)
    # wrap a function that can build the ODE problems
    function _problem_builder(x, v)
        geodesic_ode_problem(trace, config.metric, x, v, config.λ_domain, cbs)
    end
    wrap_arguments(config, _problem_builder, trace.μ)
end

wrap_arguments(config, f, μ) =
    wrap_arguments(config.metric, config.position, config.velocity, f, μ)

function wrap_arguments(
    m::AbstractMetric,
    init_pos::U,
    init_vel::V,
    _problem_func,
    μ,
) where {U,V}
    if U <: SVector && V <: SVector
        # single position and velocity
        return _problem_func(init_pos, constrain_all(m, init_pos, init_vel, μ))
    elseif U <: SVector && V <: Function
        # single position, velocity generating function
        _vfunc = wrap_constraint(m, init_pos, init_vel, μ)
        prob = _problem_func(init_pos, _vfunc(1))
        ens_prob = EnsembleProblem(
            prob,
            prob_func = (prob, i, repeat) -> _problem_func(init_pos, _vfunc(i)),
            safetycopy = false,
        )
        return ens_prob
    elseif U === V
        # both are arrays of SVectors
        _vels = constrain_all(m, init_pos, init_vel, μ)
        prob = _problem_func(init_pos[1], _vels[1])
        ens_prob = EnsembleProblem(
            prob,
            prob_func = (prob, i, repeat) -> _problem_func(init_pos[i], _vels[i]),
            safetycopy = false,
        )
        return ens_prob
    else
        error("No method for wrapping position and velocity types: $U, $V")
    end
end
