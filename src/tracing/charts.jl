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
            integrator.p.status = StatusCodes.WithinInnerBoundary
            terminate!(integrator)
        else
            integrator.p.status = StatusCodes.OutOfDomain
            terminate!(integrator)
        end
    end
    DiscreteCallback(on_chart, chart_terminate!)
end

function chart_for_metric(
    m::AbstractMetricParams{T},
    outer_radius = 1200.0;
    closest_approach = 1.01,
) where {T}
    chart =
        PolarChart(GradusBase.inner_radius(m) * closest_approach, convert(T, outer_radius))
    chart
end
