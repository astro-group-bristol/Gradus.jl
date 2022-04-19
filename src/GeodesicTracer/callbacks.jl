function metric_callback(m::AbstractMetricParams{T}) where {T}
    min_r = inner_radius(m)
    DiscreteCallback((u, λ, integrator) -> u[6] ≤ min_r * 1.1 || u[6] > 1200.0, terminate!)
end

function create_callback_set(m::AbstractMetricParams{T}, cb::Nothing) where {T}
    mcb = metric_callback(m)
    mcb isa Tuple ? CallbackSet(mcb...) : mcb
end

function create_callback_set(
    m::AbstractMetricParams{T},
    cb::NTuple{N,SciMLBase.DECallback},
) where {T,N}
    mcb = metric_callback(m)
    mcb isa Tuple ? CallbackSet(mcb..., cb...) : CallbackSet(mcb, cb...)
end

function create_callback_set(m::AbstractMetricParams{T}, cb::SciMLBase.DECallback) where {T}
    mcb = metric_callback(m)
    mcb isa Tuple ? CallbackSet(mcb..., cb) : CallbackSet(mcb, cb)
end
