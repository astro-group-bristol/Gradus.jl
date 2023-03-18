@inline function tracing_configuration(
    m::AbstractMetricParameters,
    position,
    velocity,
    geometry::AbstractAccretionGeometry,
    args...;
    gtol = 1e-2,
    callback = nothing,
    kwargs...,
)
    _tracing_configuration(
        m,
        position,
        velocity,
        geometry,
        args...;
        callback = add_collision_callback(callback, geometry; gtol = gtol),
        kwargs...,
    )
end

function add_collision_callback(callback, accretion_geometry; gtol)
    geometry_cb = build_collision_callback(accretion_geometry; gtol = gtol)
    merge_callbacks(callback, geometry_cb)
end

"""
    build_collision_callback(m::AbstractAccretionGeometry{T})

Generates the callback used for the integration. Returns a `Function`, with the fingerprint
```julia
function callback(u, λ, integrator)::Bool
    # ...
end
```
"""
function build_collision_callback(g::AbstractAccretionGeometry; gtol)
    DiscreteCallback(
        (u, λ, integrator) ->
            intersects_geometry(g, line_element(u, integrator), integrator),
        terminate_with_status!(StatusCodes.IntersectedWithGeometry),
    )
end
