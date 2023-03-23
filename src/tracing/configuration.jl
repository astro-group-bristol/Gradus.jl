
struct TracingConfiguration{
    T,
    MetricType,
    PositionType,
    VelocityType,
    GeometryType,
    CallbackType,
    ChartType,
    SolverType,
    EnsembleType,
    TrajectoriesType,
}
    metric::MetricType
    position::PositionType
    velocity::VelocityType
    geometry::GeometryType
    chart::ChartType
    callback::CallbackType
    # numerics
    solver::SolverType
    ensemble::EnsembleType
    trajectories::TrajectoriesType
    λ_domain::Tuple{T,T}
    abstol::T
    reltol::T
    # constructor
    function TracingConfiguration(
        m::AbstractMetric{T},
        position,
        velocity::V,
        geometry,
        chart,
        callback,
        solver,
        ensemble,
        trajectories,
        λ_min,
        λ_max,
        abstol,
        reltol,
    ) where {T,V}
        if V <: Function && isnothing(trajectories)
            error("When velocity is a function, trajectories must be defined.")
        end
        if V <: SVector && !isnothing(trajectories)
            error(
                "Trajectories should be `nothing` when solving only a single geodesic problem.",
            )
        end
        _trajectories = eltype(velocity) <: SVector ? length(velocity) : trajectories
        _ensemble = restrict_ensemble(m, ensemble)
        new{
            T,
            typeof(m),
            typeof(position),
            typeof(velocity),
            typeof(geometry),
            typeof(callback),
            typeof(chart),
            typeof(solver),
            typeof(_ensemble),
            typeof(_trajectories),
        }(
            m,
            position,
            velocity,
            geometry,
            chart,
            callback,
            solver,
            _ensemble,
            _trajectories,
            (λ_min, λ_max),
            abstol,
            reltol,
        )
    end
end

@inline function _tracing_configuration(
    m::AbstractMetric,
    position,
    velocity,
    geometry,
    λs::D,
    ;
    chart = chart_for_metric(m),
    callback = nothing,
    solver = Tsit5(),
    ensemble = EnsembleThreads(),
    trajectories = nothing,
    abstol = 1e-9,
    reltol = 1e-9,
    solver_opts...,
) where {D}
    _velocity, _trajectories = promote_velfunc(m, position, velocity, trajectories)
    config = TracingConfiguration(
        m,
        position,
        _velocity,
        geometry,
        chart,
        callback,
        solver,
        ensemble,
        _trajectories,
        D <: Number ? D(0) : λs[1],
        D <: Number ? λs : λs[2],
        abstol,
        reltol,
    )
    config, solver_opts
end

@inline function tracing_configuration(
    ::AbstractTraceParameters,
    m::AbstractMetric,
    position,
    velocity,
    λs;
    kwargs...,
)
    _tracing_configuration(m, position, velocity, nothing, λs; kwargs...)
end

promote_velfunc(m, position, velocity, trajectories) = (velocity, trajectories)
