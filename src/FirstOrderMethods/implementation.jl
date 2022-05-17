
@with_kw struct FirstOrderGeodesicPoint{T,V,P} <: AbstractGeodesicPoint{T}
    retcode::Symbol
    "Start time"
    t1::T
    "End time"
    t2::T
    "Start position"
    u1::V
    "End position"
    u2::V
    "Start velocity"
    v1::V
    "End velocity"
    v2::V
    "First order extra parameter"
    p::P
end

@inbounds function getgeodesicpoint(m::AbstractFirstOrderMetricParams{T}, sol::SciMLBase.AbstractODESolution{T,N,S}) where {T,N,S}
    us, ts, p = unpack_solution(sol)

    u_start = SVector{4,T}(us[1][1:4])
    v_start = SVector{4,T}(us[1][5:8])
    t_start = SVector{4,T}(ts[1])

    u_end = SVector{4,T}(us[end][1:4])
    v_end = SVector{4,T}(us[end][5:8])
    t_end = SVector{4,T}(ts[end])

    FirstOrderGeodesicPoint(sol.retcode, t_start, t_end, u_start, u_end, v_start, v_end, p)
end

function metric_callback(
    m::AbstractFirstOrderMetricParams{T},
    closest_approach,
    effective_infinity,
) where {T}
    (
        ensure_chart_callback(m, closest_approach, effective_infinity),
        DiscreteCallback(radial_negative_check(m), flip_radial_sign!),
        DiscreteCallback(angular_negative_check(m), flip_angular_sign!),
    )
end

four_velocity(u, m::AbstractFirstOrderMetricParams{T}, p) where {T} =
    error("Not implmented for $(typeof(m)).")

make_parameters(L, Q, sign_θ, T) =
    (L = L, Q = Q, r = -1, θ = convert(Int, sign_θ), changes = T[0.0, 0.0])

function integrator_problem(
    m::AbstractFirstOrderMetricParams{T},
    pos::StaticVector{S,T},
    vel::StaticVector{S,T},
    time_domain,
) where {S,T}
    L, Q = calc_lq(m, pos, vel)
    ODEProblem{false}(pos, time_domain, make_parameters(L, Q, vel[2], T)) do u, p, λ
        SVector(four_velocity(u, m, p)...)
    end
end

function integrator_problem(
    m::AbstractFirstOrderMetricParams{T},
    pos::AbstractVector{T},
    vel::AbstractVector{T},
    time_domain,
) where {T}
    L, Q = calc_lq(m, pos, vel)
    ODEProblem{true}(pos, time_domain, make_parameters(L, Q, vel[2], T)) do du, u, p, λ
        du .= four_velocity(u, m, p)
    end
end

convert_velocity_type(u::StaticVector{S,T}, v) where {S,T} = convert(SVector{S,T}, v)
convert_velocity_type(u::AbstractVector{T}, v) where {T} = convert(typeof(u), collect(v))

Vr(m::AbstractFirstOrderMetricParams{T}, u, p) where {T} =
    error("Not implmented for $(typeof(m)).")
Vθ(m::AbstractFirstOrderMetricParams{T}, u, p) where {T} =
    error("Not implmented for $(typeof(m)).")

calc_lq(m::AbstractFirstOrderMetricParams{T}, pos, vel) where {T} =
    error("Not implmented for $(typeof(m)).")
