@inline function tracing_configuration(
    trace::AbstractTrace,
    m::AbstractMetric,
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
function geometry_collision_callback(g::AbstractAccretionGeometry, ::AbstractTrace; gtol)
    DiscreteCallback(
        (u, λ, integrator) ->
            intersects_geometry(g, line_element(u, integrator), integrator),
        terminate_with_status!(StatusCodes.IntersectedWithGeometry),
    )
end

function geometry_collision_callback(
    g::AbstractAccretionDisc{T},
    ::AbstractTrace;
    gtol,
    interp_points = 12,
) where {T}
    ContinuousCallback(
        (u, λ, integrator) -> distance_to_disc(g, u; gtol = gtol),
        terminate_with_status!(StatusCodes.IntersectedWithGeometry),
        interp_points = interp_points,
        save_positions = (false, false),
    )
end

function geometry_collision_callback(
    cg::CompositeGeometry{T},
    ::AbstractTraceParameters;
    gtol,
    interp_points = 8,
) where {T}
    map(cg.geometry) do g
        ContinuousCallback(
            (u, λ, integrator) -> distance_to_disc(g, u; gtol = gtol),
            terminate_with_status!(StatusCodes.IntersectedWithGeometry),
            interp_points = interp_points,
            save_positions = (true, false),
        )
    end
end

function intersected_with_geometry(gps::AbstractArray{<:AbstractGeodesicPoint}, limiter)
    [(i.status == StatusCodes.IntersectedWithGeometry) && limiter(i.x) for i in gps]
end
