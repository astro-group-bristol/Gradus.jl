struct IntegrationParameters{M} <: AbstractIntegrationParameters{M}
    metric::M
    status::MutStatusCode
    IntegrationParameters(metric::M, status) where {M} =
        new{M}(metric, MutStatusCode(status))
end

set_status_code!(params::IntegrationParameters, status::StatusCodes.T) =
    params.status[1] = status
get_status_code(params::IntegrationParameters) = params.status[1]
get_metric(params::IntegrationParameters) = params.metric

function Base.show(io::IO, @nospecialize(p::AbstractIntegrationParameters))
    status = get_status_code(p)
    name = Base.typename(typeof(p)).name
    print(io, "$name(status=$status)")
end

"""
    geodesic_ode_problem(
        trace::AbstractTrace,
        m::AbstractMetric,
        pos, 
        vel,
        time_domain::Tuple,
        callback
    

Returns an `OrdinaryDiffEq.ODEProblem{false}`, specifying the ODE problem to be solved. 
The precise problem depends on the [`AbstractTrace`](@ref) and [`AbstractMetric`](@ref) defined.

May be overwritten to more easily define a new tracing problem. The standard geodesic equation implemention looks like:

```julia
function geodesic_ode_problem(
    ::TraceGeodesic,
    m::AbstractMetric,
    pos,
    vel,
    time_domain,
    callback,
)
    function f(u::SVector{8,T}, p, λ) where {T}
        @inbounds let x = SVector{4,T}(@view(u[1:4])), 
            v = SVector{4,T}(@view(u[5:8]))
            dv = SVector{4,T}(geodesic_equation(m, x, v))
            # modify the differential equation here
            vcat(v, dv)
        end
    end
    
    # add additional parameters here
    u_init = vcat(pos, vel)
    ODEProblem{false}(
        f,
        u_init,
        time_domain,
        # specify parameters needed by `f` here
        IntegrationParameters(StatusCodes.NoStatus);
        callback = callback,
    )
end
```

See also [`TraceGeodesic`](@ref) and [`TraceRadiativeTransfer`](@ref).
"""
function geodesic_ode_problem(
    ::TraceGeodesic,
    m::AbstractMetric,
    pos,
    vel,
    time_domain,
    callback,
)
    u_init = vcat(pos, vel)
    ODEProblem{false}(
        _second_order_ode_f,
        u_init,
        time_domain,
        IntegrationParameters(m, StatusCodes.NoStatus);
        callback = callback,
    )
end

function _second_order_ode_f(u::SVector{8}, p, λ)
    x = @inbounds SVector{4}(view(u, 1:4))
    v = @inbounds SVector{4}(view(u, 5:8))
    dv = geodesic_equation(get_metric(p), x, v)
    return vcat(v, dv)
end

"""
    assemble_tracing_problem(trace::AbstractTrace, config::TracingConfiguration)

Merges callbacks, defines an ODE builder through (a variation of) [`geodesic_ode_problem`](@ref),
and wraps the ODE problem depending on the input argument types. 

This function need only be overwritten if the [`AbstractTrace`](@ref) requires fine control
or non-standard arguments when building the ODE. See, e.g., the [`TraceRadiativeTransfer`](@ref)
implementation.

For merging the callbacks, use [`create_callback_set`](@ref).

For wrapping arguments, use the utility function [`wrap_arguments`](@ref).
"""
function assemble_tracing_problem(trace::AbstractTrace, config::TracingConfiguration)
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
