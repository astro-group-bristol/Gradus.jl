@inline function tracing_configuration(
    trace::AbstractTrace,
    m::AbstractMetric,
    position,
    velocity,
    geometry::AbstractAccretionGeometry,
    args...;
    gtol = 1e-4,
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
        save_positions = (false, false),
    )
end

function geometry_collision_callback(
    g::AbstractAccretionDisc,
    trace::AbstractTrace;
    gtol,
    interp_points = 8,
)
    ContinuousCallback(
        intersection_callbacks(g, trace, 1; gtol = gtol)...,
        interp_points = interp_points,
        save_positions = (false, false),
    )
end

function _intersection_condition(g::AbstractAccretionDisc; gtol)
    function _distance_to_disc_wrapper(u, λ, integrator)
        distance_to_disc(g, u; gtol = gtol)
    end
end

function _intersection_affect(::AbstractAccretionDisc)
    terminate_with_status!(StatusCodes.IntersectedWithGeometry)
end

"""
    intersection_callbacks

Used to control the trace dispatch.
"""
function intersection_callbacks(g::AbstractAccretionGeometry, ::AbstractTrace, ::Int; gtol)
    _intersection_condition(g; gtol = gtol), _intersection_affect(g)
end

"""
    geometry_collision_callback(cg::CompositeGeometry, trace::AbstractTrace; kwargs...)

Assemble a `VectorContinuousCallback`.
"""
function geometry_collision_callback(
    cg::CompositeGeometry,
    trace::AbstractTrace;
    gtol,
    interp_points = 8,
)
    # geometry counter so we can associate a unique id with each piece of geometry
    i = 0

    callbacks = map(cg.geometry) do g
        i += 1
        intersection_callbacks(g, trace, i; gtol = gtol)
    end

    function _composite_condition(out, u, t, integ)
        results = map(callbacks) do f
            f[1](u, t, integ)
        end

        for i = 1:length(results)
            out[i] = results[i]
        end
    end

    function _composite_affect!(integ, idx)
        callbacks[idx][2](integ)
    end

    VectorContinuousCallback(
        _composite_condition,
        _composite_affect!,
        length(cg.geometry),
        interp_points = interp_points,
    )
end

function intersected_with_geometry(gps::AbstractArray{<:AbstractGeodesicPoint}, limiter)
    [(i.status == StatusCodes.IntersectedWithGeometry) && limiter(i.x) for i in gps]
end
