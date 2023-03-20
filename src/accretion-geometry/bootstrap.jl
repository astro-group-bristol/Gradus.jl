@inline function tracing_configuration(
    trace::AbstractTraceParameters,
    m::AbstractMetricParameters,
    position,
    velocity,
    geometry::AbstractAccretionGeometry,
    args...;
    gtol = 1e-2,
    callback = nothing,
    kwargs...,
)
    geometry_callback = geometry_collision_callback(geometry, trace, gtol = gtol)
    _tracing_configuration(
        m,
        position,
        velocity,
        geometry,
        args...;
        callback = merge_callbacks(callback, geometry_callback),
        kwargs...,
    )
end

"""
    geometry_collision_callback(m::AbstractAccretionGeometry{T})

Generates the callback used for the integration. Returns a `Function`, with the fingerprint
```julia
function callback(u, λ, integrator)::Bool
    # ...
end
```
"""
function geometry_collision_callback(
    g::AbstractAccretionGeometry,
    ::AbstractTraceParameters;
    gtol,
)
    DiscreteCallback(
        (u, λ, integrator) ->
            intersects_geometry(g, line_element(u, integrator), integrator),
        terminate_with_status!(StatusCodes.IntersectedWithGeometry),
    )
end

function geometry_collision_callback(
    g::AbstractAccretionDisc{T},
    ::AbstractTraceParameters;
    gtol,
    interp_points = 8,
) where {T}
    ContinuousCallback(
        (u, λ, integrator) -> distance_to_disc(g, u; gtol = gtol),
        terminate_with_status!(StatusCodes.IntersectedWithGeometry),
        interp_points = interp_points,
        save_positions = (true, false),
    )
end


function intersected_with_geometry(gps::AbstractArray{<:AbstractGeodesicPoint}, limiter)
    [(i.status == StatusCodes.IntersectedWithGeometry) && limiter(i.x) for i in gps]
end
