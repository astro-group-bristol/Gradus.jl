function terminate_with_status!(status::StatusCodes.T)
    integrator -> begin
        integrator.p.status = status
        terminate!(integrator)
    end
end

@inline function ensure_chart_callback(
    m::AbstractMetricParams,
    closest_approach,
    effective_infinity,
)
    min_r = inner_radius(m)
    # terminate integration if we come within some % of the black hole radius
    DiscreteCallback(
        (u, λ, integrator) -> u[2] ≤ min_r * closest_approach || u[2] > effective_infinity,
        terminate_with_status!(StatusCodes.OutOfDomain),
    )
end

function metric_callback(m::AbstractMetricParams, closest_approach, effective_infinity)
    ensure_chart_callback(m, closest_approach, effective_infinity)
end

function create_callback_set(
    m::AbstractMetricParams,
    cb::C,
    closest_approach,
    effective_infinity,
) where {C}
    mcb = metric_callback(m, closest_approach, effective_infinity)
    if C <: SciMLBase.DECallback
        mcb isa Tuple ? CallbackSet(cb, mcb...) : CallbackSet(cb, mcb)
    elseif C <: Union{Tuple,AbstractVector}
        mcb isa Tuple ? CallbackSet(cb..., mcb...) : CallbackSet(cb..., mcb)
    else
        # else C is Nothing
        mcb isa Tuple ? CallbackSet(mcb...) : mcb
    end
end

# predefined callbacks
function domain_upper_hemisphere(δ = 1e-3)
    DiscreteCallback((u, t, integrator) -> u[2] * cos(u[3]) < δ, terminate_with_status!(StatusCodes.OutOfDomain))
end

export domain_upper_hemisphere
