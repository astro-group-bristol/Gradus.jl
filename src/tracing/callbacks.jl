function terminate_with_status!(status::StatusCodes.T)
    function _terminate_with_status_closure!(integrator)
        integrator.p.status = status
        terminate!(integrator)
    end
end

# this function gets specialised for different metric parameters types, e.g. for first order
@inline function metric_callback(::AbstractMetricParams, chart)
    chart_callback(chart)
end

function create_callback_set(m::AbstractMetricParams, cb::C, chart) where {C}
    mcb = metric_callback(m, chart)
    if C <: SciMLBase.DECallback
        mcb isa Tuple ? CallbackSet(cb, mcb...) : CallbackSet(cb, mcb)
    elseif C <: Tuple
        mcb isa Tuple ? CallbackSet(cb..., mcb...) : CallbackSet(cb..., mcb)
    elseif C <: Nothing
        mcb isa Tuple ? CallbackSet(mcb...) : mcb
    else
        error("Unknown callback type $C")
    end
end

# predefined callbacks
function domain_upper_hemisphere(δ = 1e-3)
    DiscreteCallback(
        (u, t, integrator) -> u[2] * cos(u[3]) < δ,
        terminate_with_status!(StatusCodes.OutOfDomain),
    )
end

export domain_upper_hemisphere
