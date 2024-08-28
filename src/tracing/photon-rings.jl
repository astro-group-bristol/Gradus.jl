function winding_callback(inc)
    function _domain_upper_hemisphere_check(u, t, integrator)
        if (integrator.p.winding[] % 2 == 0)
            # even number of windings
            u[3] > inc
        else
            # odd number
            u[3] < inc
        end
    end
    function _bump_winding_number!(integrator)
        integrator.p.winding[] += 1
    end
    DiscreteCallback(_domain_upper_hemisphere_check, _bump_winding_number!)
end

struct TraceWindings{T} <: AbstractTrace
    "Geodesic mass"
    μ::T
    "Plane through which to count the intersections. Defaults to equitorial (π/2)."
    plane_inc::T
end
TraceWindings() = TraceWindings(0.0, π / 2)

const MutWinding = MVector{1,Int}
struct WindingParams{M} <: AbstractIntegrationParameters{M}
    metric::M
    status::MutStatusCode
    winding::MutWinding
    WindingParams(metric::M, status) where {M} =
        new{M}(metric, MutStatusCode(status), MutWinding(0))
end

function update_integration_parameters!(p::WindingParams, new::WindingParams)
    set_status_code!(p, get_status_code(new))
    p.winding .= new.winding
    p
end

function merge_auxiliary(p::WindingParams, aux)
    if isnothing(aux)
        (; winding = p.winding[])
    else
        (; winding = p.winding[], aux)
    end
end

set_status_code!(params::WindingParams, status::StatusCodes.T) = params.status[1] = status
get_status_code(params::WindingParams) = params.status[1]
get_metric(params::WindingParams) = params.metric

function geodesic_ode_problem(
    trace::TraceWindings,
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
        WindingParams(m, StatusCodes.NoStatus);
        callback = merge_callbacks(callback, winding_callback(trace.plane_inc)),
    )
end


export TraceWindings
