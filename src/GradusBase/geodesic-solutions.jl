abstract type AbstractGeodesicPoint{T} end

@with_kw struct GeodesicPoint{T} <: AbstractGeodesicPoint{T}
    retcode::Symbol
    t::T
    u::AbstractVector{T}
    v::AbstractVector{T}
    # we don't store the problem parameters
    # and can create a specialistion for the carter method
    # then provide dispatched accessor methods
    # p::P
end

function geodesic_point_type(m::AbstractMetricParams{T}) where {T}
    GeodesicPoint{T}
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

function get_endpoint(
    m::AbstractMetricParams{T},
    sol::SciMLBase.AbstractODESolution{T,N,S},
) where {T,N,S}
    us, ts, _ = unpack_solution(sol)
    u = us[end][1:4]
    v = us[end][5:8]
    t = ts[end]
    GeodesicPoint(sol.retcode, t, u, v)
end

function get_endpoint(
    m::AbstractMetricParams{T},
    simsol::SciMLBase.AbstractEnsembleSolution{T,N,S},
) where {T,N,S}
    map(sol -> get_endpoint(m, sol), simsol)
end
