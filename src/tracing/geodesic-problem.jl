function geodesic_ode_problem(
    m::AbstractMetricParams,
    pos::StaticVector,
    vel::StaticVector,
    time_domain;
    q = 0.0,
    kwargs...,
)
    function f(u::SVector{8,T}, p, λ) where {T}
        @inbounds let x = SVector{4,T}(@view(u[1:4])), v = SVector{4,T}(@view(u[5:8]))
            dv = SVector{4,T}(geodesic_eq(m, x, v))
            # SVector{8}(v[1], v[2], v[3], v[4], dv[1], dv[2], dv[3], dv[4])
            vcat(v, dv)
        end
    end

    u_init = vcat(pos, vel)
    ODEProblem{false}(
        f,
        u_init,
        time_domain,
        IntegrationParameters(StatusCodes.NoStatus);
        kwargs...,
    )
end

function geodesic_problem(
    m::AbstractMetricParams,
    init_pos::U,
    init_vel::V,
    time_domain::NTuple{2},
    ;
    callback = nothing,
    chart = chart_for_metric(m),
    μ = 0.0,
    kwargs...,
) where {U,V}
    # create the callback set for the problem
    cbs = create_callback_set(m, callback, chart)
    function _problem(x, v)
        geodesic_ode_problem(m, x, v, time_domain; callback = cbs, kwargs...)
    end
    wrap_arguments(m, init_pos, init_vel, _problem, μ)
end

function wrap_arguments(
    m::AbstractMetricParams,
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
    end
end
