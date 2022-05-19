abstract type AbstractGeodesicPoint{T} end

@with_kw struct GeodesicPoint{T,V<:AbstractVector} <: AbstractGeodesicPoint{T}
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
    # we don't store the problem parameters
    # and can create a specialistion for the carter method
    # then provide dispatched accessor methods
    # p::P
end

@inbounds function getgeodesicpoint(
    m::AbstractMetricParams{T},
    sol::SciMLBase.AbstractODESolution{T,N,S},
) where {T,N,S}
    us, ts, _ = unpack_solution(sol)

    u_start = SVector{4,T}(us[1][1:4])
    v_start = SVector{4,T}(us[1][5:8])
    t_start = ts[1]

    u_end = SVector{4,T}(us[end][1:4])
    v_end = SVector{4,T}(us[end][5:8])
    t_end = ts[end]

    GeodesicPoint(sol.retcode, t_start, t_end, u_start, u_end, v_start, v_end)
end

# TODO: GeodesicPath structure for the full geodesic path
# do we want to support this?

function unpack_solution(sol::SciMLBase.AbstractODESolution{T,N,S}) where {T,N,S}
    u = sol.u
    p = sol.prob.p
    t = sol.t
    (u, t, p)
end

function unpack_solution(simsol::SciMLBase.AbstractEnsembleSolution{T,N,V}) where {T,N,V}
    map(unpack_solution, simsol)
end
