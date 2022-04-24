
@with_kw struct FirstOrderGeodesicPoint{T,P} <: AbstractGeodesicPoint{T}
    retcode::Symbol
    t::T
    u::AbstractVector{T}
    v::AbstractVector{T}
    p::P
end

function geodesic_point_type(m::AbstractFirstOrderMetricParams{T}) where {T}
    p_type = typeof(make_parameters(T(0), T(0), 1, T))
    FirstOrderGeodesicPoint{T,p_type}
end

function get_endpoint(
    m::AbstractFirstOrderMetricParams{T},
    sol::SciMLBase.AbstractODESolution{T,N,S},
) where {T,N,S}
    us, ts, p = unpack_solution(sol)
    u = us[end]
    v = four_velocity(u, m, p)
    t = ts[end]
    FirstOrderGeodesicPoint(sol.retcode, t, u, convert_velocity_type(u, v), p)
end

function metric_callback(
    m::AbstractFirstOrderMetricParams{T};
    closest_approach = 1.01,
    effective_infinity = 1200.0,
) where {T}
    (
        DiscreteCallback(
            ensure_domain(m, closest_approach, effective_infinity),
            terminate!,
        ),
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
