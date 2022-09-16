"""
    AbstractFirstOrderMetricParams{T} <: AbstractMetricParams{T}

Abstract type for metrics using the 1st-order integration method. The 1st-order methods
reuse the velocity vector as a parameter vector, where only element `vel[2]` and `vel[3]`
are used, and are local observer ratios ``\\sin \\Theta`` and ``\\sin \\Phi`` respectively.

Require implementation of
- [`inner_radius`](@ref)
- [`constrain`](@ref)
- [`four_velocity`](@ref)
- [`calc_lq`](@ref)
- [`Vr`](@ref)
- [`Vθ`](@ref)
- [`alpha_beta_to_vel`](@ref)
"""
abstract type AbstractFirstOrderMetricParams{T} <: AbstractMetricParams{T} end

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

@inbounds function getgeodesicpoint(
    m::AbstractFirstOrderMetricParams{T},
    sol::SciMLBase.AbstractODESolution{T,N,S},
) where {T,N,S}
    us, ts, p = unpack_solution(sol)

    u_start = SVector{4,T}(us[1][1:4])
    v_start = SVector{4,T}(four_velocity(u_start, m, p))
    t_start = ts[1]

    u_end = SVector{4,T}(us[end][1:4])
    v_end = SVector{4,T}(four_velocity(u_end, m, p))
    t_end = ts[end]

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

"""
    $(TYPEDSIGNATURES)

Calculate the four-velocity at a point `u`, given a set of metric parameters and the constants
of motion in `p`.
"""
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

"""
    $(TYPEDSIGNATURES)

Effective potential in the radial direction. Used only to track sign changes.
"""
Vr(m::AbstractFirstOrderMetricParams{T}, u, p) where {T} =
    error("Not implmented for $(typeof(m)).")
"""
    $(TYPEDSIGNATURES)

Effective potential in the angular direction. Used only to track sign changes.
"""
Vθ(m::AbstractFirstOrderMetricParams{T}, u, p) where {T} =
    error("Not implmented for $(typeof(m)).")

"""
    $(TYPEDSIGNATURES)

Calculate constants of motion ``L`` and ``Q``, given a set of metric parameters,
the geodesic position, and the `param` vector.
"""
calc_lq(m::AbstractFirstOrderMetricParams{T}, pos, param) where {T} =
    error("Not implmented for $(typeof(m)).")

function flip_radial_sign!(integrator)
    p = integrator.p
    integrator.p = @set p.r = -p.r
    integrator.sol.prob.p.changes[1] = integrator.t[end]
end

function flip_angular_sign!(integrator)
    p = integrator.p
    integrator.p = @set p.θ = -p.θ
    integrator.sol.prob.p.changes[2] = integrator.t[end]
end

function radial_negative_check(m::AbstractFirstOrderMetricParams{T}) where {T}
    (u, λ, integrator) -> Vr(m, u, integrator.p) < 0
end

function angular_negative_check(m::AbstractFirstOrderMetricParams{T}) where {T}
    (u, λ, integrator) -> Vθ(m, u, integrator.p) < 0
end
