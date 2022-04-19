

@with_kw struct FirstOrderGeodesicPoint{T,P} <: AbstractGeodesicPoint{T}
    retcode::Symbol
    t::T
    u::AbstractVector{T}
    v::AbstractVector{T}
    p::P
end

function geodesic_point_type(m::CarterMethodBL{T}) where {T}
    p_type = typeof(make_parameters(T(0), T(0), 1, T))
    FirstOrderGeodesicPoint{T,p_type}
end

make_parameters(L, Q, sign_θ, T) =
    (L = L, Q = Q, r = -1, θ = convert(Int, sign_θ), changes = T[0.0, 0.0])


function metric_callback(m::CarterMethodBL{T}) where {T}
    (
        DiscreteCallback(ensure_domain(m), terminate!),
        DiscreteCallback(radial_negative_check(m), flip_radial_sign!),
        DiscreteCallback(angular_negative_check(m), flip_angular_sign!),
    )
end


convert_velocity_type(u::StaticVector{S,T}, v) where {S,T} = convert(SVector{S,T}, v)
convert_velocity_type(u::AbstractVector{T}, v) where {T} = convert(typeof(u), collect(v))

function get_endpoint(
    m::CarterMethodBL{T},
    sol::SciMLBase.AbstractODESolution{T,N,S},
) where {T,N,S}
    us, ts, p = unpack_solution(sol)
    u = us[end]
    v = carter_velocity(u, m.E, m.M, m.a, p)
    t = ts[end]
    FirstOrderGeodesicPoint(sol.retcode, t, u, convert_velocity_type(u, v), p)
end

Vr(m::AbstractFirstOrderMetricParams{T}, u, p) where {T} =
    error("Not implmented for $(typeof(m)).")
Vθ(m::AbstractFirstOrderMetricParams{T}, u, p) where {T} =
    error("Not implmented for $(typeof(m)).")


export FirstOrderGeodesicPoint
