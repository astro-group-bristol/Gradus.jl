function terminate_with_status!(status::StatusCodes.T)
    function _terminate_with_status_closure!(integrator)
        set_status_code!(integrator.p, status)
        terminate!(integrator)
    end
end

# this function gets specialised for different metric parameters types, e.g. for first order
@inline function metric_callback(::AbstractMetric, chart::AbstractChart)
    chart_callback(chart)
end

function merge_callbacks(cbs1, cb::C) where {C}
    if C <: SciMLBase.DECallback
        cbs1 isa Tuple ? CallbackSet(cb, cbs1...) : CallbackSet(cb, cbs1)
    elseif C <: Tuple
        cbs1 isa Tuple ? CallbackSet(cb..., cbs1...) : CallbackSet(cb..., cbs1)
    elseif C <: Nothing
        cbs1 isa Tuple ? CallbackSet(cbs1...) : cbs1
    else
        error("Unknown callback type $C")
    end
end

function create_callback_set(m::AbstractMetric, cb, chart::AbstractChart)
    mcb = metric_callback(m, chart)
    merge_callbacks(mcb, cb)
end

# predefined callbacks
function domain_upper_hemisphere(δ = 1e-4)
    function _domain_upper_hemisphere_check(u, t, integrator)
        _spinaxis_project(u, signed = true) < δ
    end
    DiscreteCallback(
        _domain_upper_hemisphere_check,
        terminate_with_status!(StatusCodes.OutOfDomain),
    )
end

export domain_upper_hemisphere
