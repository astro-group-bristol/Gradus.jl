"""
    AbstractFirstOrderMetric{T} <: AbstractMetric{T}

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
- [`impact_parameters_to_three_velocity`](@ref)
"""
abstract type AbstractFirstOrderMetric{T} <: AbstractMetric{T,BoyerLindquist{(:r, :θ)}} end

function restrict_ensemble(::AbstractFirstOrderMetric, ensemble::EnsembleEndpointThreads)
    @warn "First order methods do not support `EnsembleEndpointThreads`. Automatically switching to `EnsembleThreads`"
    EnsembleThreads()
end

@with_kw struct FirstOrderGeodesicPoint{T,V,P} <: AbstractGeodesicPoint{T}
    "Return code of the integrator for this geodesic."
    status::StatusCodes.T
    "Start affine time"
    λ_min::T
    "End affine time"
    λ_max::T
    "Start four position"
    x_init::V
    "End four position"
    x::V
    "Start four velocity"
    v_init::V
    "End four velocity"
    v::V
    "First order extra parameter"
    p::P
end

@inbounds function unpack_solution(
    m::AbstractFirstOrderMetric{T},
    sol::SciMLBase.AbstractODESolution{T},
) where {T}
    us, ts, p = unpack_solution(sol)

    @inbounds @views begin
        u_start = SVector{4,T}(us[1][1:4])
        v_start = SVector{4,T}(four_velocity(u_start, m, p))
        t_start = ts[1]

        u_end = SVector{4,T}(us[end][1:4])
        v_end = SVector{4,T}(four_velocity(u_end, m, p))
        t_end = ts[end]
    end

    FirstOrderGeodesicPoint(
        p.status,
        t_start,
        t_end,
        u_start,
        u_end,
        v_start,
        v_end,
        deepcopy(p),
    )
end

function metric_callback(m::AbstractFirstOrderMetric, chart::AbstractChart)
    (
        chart_callback(chart),
        DiscreteCallback(radial_negative_check(m), flip_radial_sign!),
        DiscreteCallback(angular_negative_check(m), flip_angular_sign!),
    )
end

"""
    $(TYPEDSIGNATURES)

Calculate the four-velocity at a point `u`, given a set of metric parameters and the constants
of motion in `p`.
"""
four_velocity(u, m::AbstractFirstOrderMetric, p) = error("Not implmented for $(typeof(m)).")

mutable struct FirstOrderIntegrationParameters{T} <: AbstractIntegrationParameters
    L::T
    Q::T
    r::Int
    θ::Int
    changes::Vector{T}
    status::StatusCodes.T
end

make_parameters(L, Q, sign_θ, ::Type{T}) where {T} =
    FirstOrderIntegrationParameters{T}(L, Q, -1, sign_θ, [0.0, 0.0], StatusCodes.NoStatus)

function update_integration_parameters!(
    p::FirstOrderIntegrationParameters,
    new::FirstOrderIntegrationParameters,
)
    p.L = new.L
    p.Q = new.Q
    p.θ = new.θ
    p.changes = new.changes
    p.status = new.status
    p
end

function geodesic_ode_problem(
    ::TraceGeodesic,
    m::AbstractFirstOrderMetric{T},
    pos::StaticVector{S,T},
    vel::StaticVector{S,T},
    time_domain,
    callback,
) where {S,T}
    L, Q = calc_lq(m, pos, vel)
    ODEProblem{false}(
        pos,
        time_domain,
        make_parameters(L, Q, vel[2], T);
        callback = callback,
    ) do u, p, λ
        SVector(four_velocity(u, m, p)...)
    end
end

convert_velocity_type(::StaticVector{S,T}, v) where {S,T} = convert(SVector{S,T}, v)
convert_velocity_type(u::AbstractVector{T}, v) where {T} = convert(typeof(u), collect(v))

"""
    $(TYPEDSIGNATURES)

Effective potential in the radial direction. Used only to track sign changes.
"""
Vr(m::AbstractFirstOrderMetric{T}, u, p) where {T} =
    error("Not implmented for $(typeof(m)).")
"""
    $(TYPEDSIGNATURES)

Effective potential in the angular direction. Used only to track sign changes.
"""
Vθ(m::AbstractFirstOrderMetric{T}, u, p) where {T} =
    error("Not implmented for $(typeof(m)).")

"""
    $(TYPEDSIGNATURES)

Calculate constants of motion ``L`` and ``Q``, given a set of metric parameters,
the geodesic position, and the `param` vector.
"""
calc_lq(m::AbstractFirstOrderMetric{T}, pos, param) where {T} =
    error("Not implmented for $(typeof(m)).")

function flip_radial_sign!(integrator)
    integrator.p.r = -integrator.p.r
    integrator.p.changes[1] = integrator.t[end]
end

function flip_angular_sign!(integrator)
    integrator.p.θ = -integrator.p.θ
    integrator.p.changes[2] = integrator.t[end]
end

function radial_negative_check(m::AbstractFirstOrderMetric{T}) where {T}
    (u, λ, integrator) -> Vr(m, u, integrator.p) < 0
end

function angular_negative_check(m::AbstractFirstOrderMetric{T}) where {T}
    (u, λ, integrator) -> Vθ(m, u, integrator.p) < 0
end

export AbstractFirstOrderMetric, FirstOrderGeodesicPoint
