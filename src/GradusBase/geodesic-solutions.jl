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

function getpoint(
    m::AbstractMetricParams{T},
    sol::SciMLBase.AbstractODESolution{T,N,S},
    i::Int,
) where {T,N,S}
    us, ts, _ = unpack_solution(sol)
    if i == -1
        u = us[end][1:4]
        v = us[end][5:8]
        t = ts[end]
        GeodesicPoint(sol.retcode, t, u, v)
    else
        u = us[i][1:4]
        v = us[i][5:8]
        t = ts[i]
        GeodesicPoint(sol.retcode, t, u, v)
    end
end

function getpoint(
    m::AbstractMetricParams{T},
    simsol::Union{SciMLBase.AbstractEnsembleSolution{T,N,S},AbstractArray{P}},
    i,
) where {T,N,S,P<:SciMLBase.AbstractODESolution}
    map(sol -> getpoint(m, sol, i), simsol)
end

function getendpoint(m::AbstractMetricParams{T}, sol_or_simsols) where {T}
    getpoint(m, sol_or_simsols, -1)
end

function getstartpoint(m::AbstractMetricParams{T}, sol_or_simsols) where {T}
    getpoint(m, sol_or_simsols, 1)
end
