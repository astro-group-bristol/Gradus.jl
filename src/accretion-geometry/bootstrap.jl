function tracegeodesics(
    m::AbstractMetricParams{T},
    init_positions,
    init_velocities,
    accretion_geometry,
    time_domain::Tuple{T,T};
    callback = nothing,
    gtol = 1e-2,
    kwargs...,
) where {T}
    cbs = add_collision_callback(callback, accretion_geometry; gtol = gtol)
    tracegeodesics(
        m,
        init_positions,
        init_velocities,
        time_domain;
        callback = cbs,
        kwargs...,
    )
end

function rendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    accretion_geometry,
    max_time::T;
    callback = nothing,
    gtol = 1e-2,
    kwargs...,
) where {T}
    cbs = add_collision_callback(callback, accretion_geometry; gtol = gtol)
    rendergeodesics(m, init_pos, max_time; callback = cbs, kwargs...)
end

function prerendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    accretion_geometry,
    max_time::T;
    callback = nothing,
    gtol = 1e-2,
    kwargs...,
) where {T}
    cbs = add_collision_callback(callback, accretion_geometry; gtol = gtol)
    prerendergeodesics(m, init_pos, max_time; callback = cbs, kwargs...)
end

add_collision_callback(::Nothing, accretion_geometry; gtol) =
    build_collision_callback(accretion_geometry; gtol = gtol)
add_collision_callback(callback::Base.AbstractVecOrTuple, accretion_geometry; gtol) =
    (callback..., build_collision_callback(accretion_geometry; gtol = gtol))
add_collision_callback(callback, accretion_geometry; gtol) =
    (callback, build_collision_callback(accretion_geometry; gtol = gtol))

"""
    build_collision_callback(m::AbstractAccretionGeometry{T})

Generates the callback used for the integration. Returns a `Function`, with the fingerprint
```julia
function callback(u, λ, integrator)::Bool
    # ...
end
```
"""
function build_collision_callback(g::AbstractAccretionGeometry{T}; gtol) where {T}
    DiscreteCallback(
        (u, λ, integrator) ->
            intersects_geometry(g, line_element(u, integrator), integrator),
        terminate_with_status!(StatusCodes.IntersectedWithGeometry),
    )
end
