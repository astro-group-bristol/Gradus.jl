function trace_single_orbit(m, r, vϕ; max_time = 300.0, μ = 1.0, θ₀ = π / 2, tracer_args...)
    # fixed in equatorial plane
    u = @SVector [0.0, r, θ₀, 0.0]
    v = @SVector [0.0, 0.0, 0.0, vϕ]
    tracegeodesics(m, u, v, (0.0, max_time); μ = μ, tracer_args...)
end

function measure_stability(m::AbstractMetric, r, vϕ; tracer_args...)
    sol = trace_single_orbit(m, r, vϕ; tracer_args...)
    rs = [sol.u[i][2] for i in eachindex(sol.u)]
    # Qs
    sum(((rs .- r) ./ r) .^ 2) / length(rs)
end

function __solve_equatorial_circular_orbit(
    m::AbstractMetric,
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

function solve_equatorial_circular_orbit(
    m::AbstractMetric,
    r::Number;
    lower_bound = 0.0,
    upper_bound = 1.0,
    optimizer = GoldenSection(),
    tracer_args...,
)
    __solve_equatorial_circular_orbit(
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

function solve_equatorial_circular_orbit(
    m::AbstractMetric,
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
        solve_equatorial_circular_orbit(
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

function trace_equatorial_circular_orbit(m::AbstractMetric, rs; kwargs...)
    map(zip(rs, solve_equatorial_circular_orbit(m, rs; kwargs...))) do (r, vϕ)
        trace_single_orbit(m, r, vϕ; kwargs...)
    end
end
function trace_equatorial_circular_orbit(m::AbstractMetric, r::Number; kwargs...)
    vϕ = solve_equatorial_circular_orbit(m, r; kwargs...)
    trace_single_orbit(m, r, vϕ; kwargs...)
end

struct PlungingInterpolation{M,_interp_type}
    m::M
    t::_interp_type
    r::_interp_type
    ϕ::_interp_type

    function PlungingInterpolation(m::M, sol) where {M<:AbstractMetric{T}} where {T}
        # sort relevant features
        I = sortperm(@view(sol[2, :]))[2:end]

        r = sol[2, :][I]

        vt = sol[5, :][I]
        vr = sol[6, :][I]
        vϕ = sol[8, :][I]

        r_interp = _make_interpolation(r, vt)
        new{M,typeof(r_interp)}(
            m,
            r_interp,
            _make_interpolation(r, vr),
            _make_interpolation(r, vϕ),
        )
    end
end

function Base.show(io::IO, ::PlungingInterpolation{M}) where {M}
    write(io, "PlungingInterpolation for $M")
end

function (pintrp::PlungingInterpolation)(r)
    r_bounded = _enforce_interpolation_bounds(r, pintrp)
    vt = pintrp.t(r_bounded)
    vr = pintrp.r(r_bounded)
    vϕ = pintrp.ϕ(r_bounded)
    SVector(vt, vr, 0, vϕ)
end

function interpolate_plunging_velocities(
    m::AbstractMetric{T};
    max_time = 50_000,
    contra_rotating = false,
    reltol = 1e-9,
    kwargs...,
) where {T}
    isco = Gradus.isco(m)

    # rule of thumb to achieve desired error
    δr = (reltol * 10)
    u = @SVector([0.0, isco - δr, deg2rad(90), 0.0])
    v = Gradus.CircularOrbits.plunging_fourvelocity(
        m,
        isco;
        contra_rotating = contra_rotating,
    )
    sol = Gradus.tracegeodesics(
        m,
        u,
        v,
        (0.0, T(max_time));
        μ = 1.0,
        reltol = reltol,
        # ensure we gets sufficiently close to the event horizon
        chart = chart_for_metric(m; closest_approach = 1.000001),
        kwargs...,
    )

    PlungingInterpolation(m, sol)
end

function _enforce_interpolation_bounds(r::Number, pintrp::PlungingInterpolation)
    _enforce_interpolation_bounds(r, first(pintrp.r.t), last(pintrp.r.t))
end

export solve_equatorial_circular_orbit,
    trace_equatorial_circular_orbit, PlungingInterpolation, interpolate_plunging_velocities
