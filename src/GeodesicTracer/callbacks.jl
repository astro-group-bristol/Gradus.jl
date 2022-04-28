@inline function ensure_chart_callback(
    m::AbstractMetricParams{T},
    closest_approach,
    effective_infinity,
) where {T}
    min_r = inner_radius(m)
    # terminate integration if we come within some % of the black hole radius
    DiscreteCallback(
        (u, λ, integrator) -> u[2] ≤ min_r * closest_approach || u[2] > effective_infinity,
        terminate!,
    )
end

function metric_callback(
    m::AbstractMetricParams{T},
    closest_approach,
    effective_infinity,
) where {T}
    ensure_chart_callback(m, closest_approach, effective_infinity)
end

function create_callback_set(
    m::AbstractMetricParams{T},
    cb::C,
    closest_approach,
    effective_infinity,
) where {T,C<:Union{NTuple{N,SciMLBase.DECallback},SciMLBase.DECallback,Nothing}} where {N}
    mcb = metric_callback(m, closest_approach, effective_infinity)
    if C <: SciMLBase.DECallback
        mcb isa Tuple ? CallbackSet(cb, mcb...) : CallbackSet(cb, mcb)
    elseif eltype(C) <: SciMLBase.DECallback
        mcb isa Tuple ? CallbackSet(cb..., mcb...) : CallbackSet(cb..., mcb)
    else
        # else C is Nothing
        mcb isa Tuple ? CallbackSet(mcb...) : mcb
    end
end
