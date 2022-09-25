"""
    abstract type AbstractGeodesicPoint

Supertype for geodesic points, used to store information about specific points along geodesic
trajectories.

!!! note
    Currently limited to storing the start and endpoint of any given trajectory. To keep the
    full geodesic path, it is encouraged to use the `SciMLBase.AbstractODESolution` directly.

Must minimally have the same fields as [`GeodesicPoint`](@ref).
Examples include [`Gradus.FirstOrderGeodesicPoint`](@ref).
"""
abstract type AbstractGeodesicPoint{T} end

"""
    struct GeodesicPoint <: AbstractGeodesicPoint

$(FIELDS)

"""
@with_kw struct GeodesicPoint{T,V<:AbstractVector} <: AbstractGeodesicPoint{T}
    "Return code of the integrator for this geodesic."
    retcode::Symbol
    "Start affine time"
    t1::T
    "End affine time"
    t2::T
    "Start four position"
    u1::V
    "End four position"
    u2::V
    "Start four velocity"
    v1::V
    "End four velocity"
    v2::V
    # we don't store the problem parameters
    # and can create a specialistion for the carter method
    # then provide dispatched accessor methods
    # p::P
end

"""
    $(TYPEDSIGNATURES)

Used to get pack a `SciMLBase.AbstractODESolution` into a [`GeodesicPoint`](@ref). Will only store
information relevant to start and endpoint of the integration. To keep the full geodesic path, it is
encouraged to use the `SciMLBase.AbstractODESolution` directly.
"""
function getgeodesicpoint(
    _::AbstractMetricParams{T},
    sol::SciMLBase.AbstractODESolution{T,N,S},
) where {T,N,S}
    @inbounds @views begin
        us, ts, _ = unpack_solution(sol)

        u_start = SVector{4,T}(us[1][1:4])
        v_start = SVector{4,T}(us[1][5:8])
        t_start = ts[1]

        u_end = SVector{4,T}(us[end][1:4])
        v_end = SVector{4,T}(us[end][5:8])
        t_end = ts[end]

        GeodesicPoint(sol.retcode, t_start, t_end, u_start, u_end, v_start, v_end)
    end
end

"""
    $(TYPEDSIGNATURES)

Unpacks each point in the solution, similar to [`getgeodesicpoint`](@ref) but returns an 
array of [`GeodesicPoint`](@ref).
"""
function getgeodesicpoints(
    _::AbstractMetricParams{T},
    sol::SciMLBase.AbstractODESolution{T,N,S},
) where {T,N,S}
    us, ts, _ = unpack_solution(sol)
    @inbounds @views begin
        u_start = SVector{4,T}(us[1][1:4])
        v_start = SVector{4,T}(us[1][5:8])
        t_start = ts[1]
        map(eachindex(us)) do i
            ui = SVector{4,T}(us[i][1:4])
            vi = SVector{4,T}(us[i][5:8])
            ti = ts[i]
            GeodesicPoint(sol.retcode, t_start, ti, u_start, ui, v_start, vi)
        end
    end
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
