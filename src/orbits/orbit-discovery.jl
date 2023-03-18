function trace_single_orbit(m, r, vϕ; max_time = 300.0, μ = 1.0, θ₀ = π / 2, tracer_args...)
    # fixed in equitorial plane
    u = @SVector [0.0, r, θ₀, 0.0]
    v = @SVector [0.0, 0.0, 0.0, vϕ]
    tracegeodesics(m, u, v, (0.0, max_time); μ = μ, tracer_args...)
end

function measure_stability(m::AbstractMetricParameters, r, vϕ; tracer_args...)
    sol = trace_single_orbit(m, r, vϕ; tracer_args...)
    rs = [sol.u[i][2] for i in eachindex(sol.u)]
    # Qs
    sum(((rs .- r) ./ r) .^ 2) / length(rs)
end

function __solve_equitorial_circular_orbit(
    m::AbstractMetricParameters,
    r,
    optimizer,
    lower_bound,
    upper_bound;
    tracer_args...,
)
    res = optimize(
        vϕ -> measure_stability(m, r, vϕ; tracer_args...),
        lower_bound,
        upper_bound,
        optimizer,
    )
    Optim.minimizer(res)
end

function solve_equitorial_circular_orbit(
    m::AbstractMetricParameters,
    r::Number;
    lower_bound = 0.0,
    upper_bound = 1.0,
    optimizer = GoldenSection(),
    tracer_args...,
)
    __solve_equitorial_circular_orbit(
        m,
        r,
        optimizer,
        lower_bound,
        upper_bound;
        tracer_args...,
    )
end

function sliding_window(func, N, lower_bound, upper_bound, lower_rate, upper_rate)
    map(1:N) do i
        res = func((i, lower_bound, upper_bound))
        lower_bound = res * lower_rate
        upper_bound = res * upper_rate
        res
    end
end

function solve_equitorial_circular_orbit(
    m::AbstractMetricParameters,
    r_range::Union{<:AbstractRange,<:AbstractArray};
    lower_bound = 0.0,
    upper_bound = 1.0,
    lower_rate = 0.88,
    upper_rate = 1.8,
    tracer_args...,
)
    r_range_reverse = sort(r_range) |> reverse
    candidate_vϕ = sliding_window(
        length(r_range_reverse),
        lower_bound,
        upper_bound,
        lower_rate,
        upper_rate,
    ) do (i, lower_bound, upper_bound)
        r = r_range_reverse[i]
        solve_equitorial_circular_orbit(
            m,
            r,
            ;
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            tracer_args...,
        )
    end
    reverse!(candidate_vϕ)
end

function trace_equitorial_circular_orbit(m::AbstractMetricParameters, rs; kwargs...)
    map(zip(rs, solve_equitorial_circular_orbit(m, rs; kwargs...))) do (r, vϕ)
        trace_single_orbit(m, r, vϕ; kwargs...)
    end
end
function trace_equitorial_circular_orbit(m::AbstractMetricParameters, r::Number; kwargs...)
    vϕ = solve_equitorial_circular_orbit(m, r; kwargs...)
    trace_single_orbit(m, r, vϕ; kwargs...)
end


export solve_equitorial_circular_orbit, trace_equitorial_circular_orbit
