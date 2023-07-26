abstract type AbstractChart end

struct PolarChart{T} <: AbstractChart
    inner_radius::T
    outer_radius::T
end

@inline function chart_callback(chart::PolarChart)
    function on_chart(u, λ, integrator)
        # terminate integration if we come within some % of the black hole radius
        u[2] ≤ chart.inner_radius || u[2] > chart.outer_radius
    end
    function chart_terminate!(integrator)
        # terminate with status code depending on whether inner or outer boundary
        if integrator.u[2] ≤ chart.inner_radius
            set_status!(integrator.p, StatusCodes.WithinInnerBoundary)
            terminate!(integrator)
        else
            set_status!(integrator.p, StatusCodes.OutOfDomain)
            terminate!(integrator)
        end
    end
    DiscreteCallback(on_chart, chart_terminate!, save_positions = (false, false))
end

struct PoloidalShapeChart{T,F} <: AbstractChart
    shapefunc::F
    outer_radius::T
end

@inline function chart_callback(chart::PoloidalShapeChart)
    function on_chart(u, λ, integrator)
        # terminate integration if we come within some % of the black hole radius
        rmin = chart.shapefunc(u[3])
        u[2] ≤ rmin || u[2] > chart.outer_radius
    end
    function chart_terminate!(integrator)
        # terminate with status code depending on whether inner or outer boundary
        rmin = chart.shapefunc(integrator.u[3])
        if integrator.u[2] ≤ rmin
            set_status!(integrator.p, StatusCodes.WithinInnerBoundary)
            terminate!(integrator)
        else
            set_status!(integrator.p, StatusCodes.OutOfDomain)
            terminate!(integrator)
        end
    end
    DiscreteCallback(on_chart, chart_terminate!)
end

function chart_for_metric(
    m::AbstractMetric{T},
    outer_radius = T(12000);
    closest_approach = T(1.01),
) where {T}
    chart =
        PolarChart(GradusBase.inner_radius(m) * closest_approach, convert(T, outer_radius))
    chart
end

function event_horizon_chart(
    m::AbstractMetric;
    outer_radius = 1200.0,
    closest_approach = 1.01,
    kwargs...,
)
    rs, θs = event_horizon(m; kwargs...)
    shapefunc = DataInterpolations.LinearInterpolation(rs .* closest_approach, θs)
    PoloidalShapeChart(shapefunc, outer_radius)
end

export event_horizon_chart
