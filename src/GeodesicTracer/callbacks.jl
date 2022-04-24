function metric_callback(m::AbstractMetricParams{T}, closest_approach, effective_infinity) where {T}
    min_r = inner_radius(m)
    # terminate integration if we come within 1% of the black hole radius
    DiscreteCallback((u, λ, integrator) -> u[6] ≤ min_r * closest_approach || u[6] > effective_infinity, terminate!)
end

function create_callback_set(m::AbstractMetricParams{T}, cb::Nothing, closest_approach, effective_infinity) where {T}
    mcb = metric_callback(m, closest_approach, effective_infinity)
    mcb isa Tuple ? CallbackSet(mcb...) : mcb
end

function create_callback_set(
    m::AbstractMetricParams{T},
    cb::NTuple{N,SciMLBase.DECallback},
    closest_approach, 
    effective_infinity
) where {T,N}
    mcb = metric_callback(m, closest_approach, effective_infinity)
    mcb isa Tuple ? CallbackSet(mcb..., cb...) : CallbackSet(mcb, cb...)
end

function create_callback_set(m::AbstractMetricParams{T}, cb::SciMLBase.DECallback, closest_approach, effective_infinity) where {T}
    mcb = metric_callback(m, closest_approach, effective_infinity)
    mcb isa Tuple ? CallbackSet(mcb..., cb) : CallbackSet(mcb, cb)
end