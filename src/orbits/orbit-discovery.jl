"""
    $(TYPEDSIGNATURES)

Quality of stability measure, which has a minima for circular orbits. Effectively
a sum of the normalised residuals.
"""
Qs(rs) = sqrt(sum((rs ./ rs[1] .- 1.0) .^ 2) / length(rs))

function trace_single_orbit(m, r, vϕ; max_time = 300.0, μ = 1.0, tracer_args...)
    # fixed in equitorial plane
    u = @SVector [0.0, r, deg2rad(90.0), 0.0]
    v = @SVector [0.0, 0.0, 0.0, vϕ]
    Gradus.tracegeodesics(m, u, v, (0.0, max_time); μ = μ, tracer_args...)
end

function measure_stability(m::AbstractMetricParams{T}, r, vϕ; tracer_args...) where {T}
    sol = trace_single_orbit(m, r, vϕ; tracer_args...)
    rs = selectdim(sol, 1, 2)
    Qs(rs)
end

function __solve_equitorial_circular_orbit(
    m::AbstractMetricParams{T},
    r,
    optimizer,
    lower_bound,
    upper_bound;
    tracer_args...,
) where {T}
    res = optimize(
        vϕ -> measure_stability(m, r, vϕ; tracer_args...),
        lower_bound,
        upper_bound,
        optimizer,
    )
    Optim.minimizer(res)
end

function solve_equitorial_circular_orbit(
    m::AbstractMetricParams{T},
    r::Number;
    lower_bound = 0.0,
    upper_bound = 1.0,
    optimizer = GoldenSection(),
    tracer_args...,
) where {T}
    __solve_equitorial_circular_orbit(m, r, optimizer, lower_bound, upper_bound)
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
    m::AbstractMetricParams{T},
    r_range::Union{AbstractRange,AbstractArray};
    lower_bound = 0.0,
    upper_bound = 1.0,
    lower_rate = 0.98,
    upper_rate = 1.5,
    kwargs...,
) where {T}
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
            lower_bound = lower_bound,
            upper_bound = upper_bound,
        )
    end
    reverse!(candidate_vϕ)
end

function trace_equitorial_circular_orbit(
    m::AbstractMetricParams{T},
    rs;
    kwargs...,
) where {T}
    map(zip(rs, solve_equitorial_circular_orbit(m, rs; kwargs...))) do (r, vϕ)
        trace_single_orbit(m, r, vϕ; kwargs...)
    end
end

export solve_equitorial_circular_orbit, trace_equitorial_circular_orbit
