abstract type AbstractIntegrationParameters end

mutable struct IntegrationParameters <: AbstractIntegrationParameters
    status::StatusCodes.T
end

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
struct GeodesicPoint{T,V<:AbstractVector} <: AbstractGeodesicPoint{T}
    "Return code of the integrator for this geodesic."
    status::StatusCodes.T
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

function Base.show(io::IO, gp::GeodesicPoint)
    print(io, "GeodesicPoint($(gp.u2...))")
end

function Base.show(io::IO, ::MIME"text/plain", gp::GeodesicPoint)
    text = """GeodesicPoint:
      . status : $(gp.status)
      . t1     : $(gp.t1)
      . t2     : $(gp.t2)
      . u1     : $(gp.u1)
      . v1     : $(gp.v1)
      . u2     : $(gp.u2)
      . v2     : $(gp.v2)
    """
    print(io, text)
end

"""
    $(TYPEDSIGNATURES)

Used to get pack a `SciMLBase.AbstractODESolution` into a [`GeodesicPoint`](@ref). Will only store
information relevant to start and endpoint of the integration. To keep the full geodesic path, it is
encouraged to use the `SciMLBase.AbstractODESolution` directly.
"""
function process_solution(_::AbstractMetricParams, sol::SciMLBase.AbstractODESolution)
    process_solution(sol)
end

function process_solution(sol::SciMLBase.AbstractODESolution{T,N,S}) where {T,N,S}
    @inbounds @views begin
        us, ts, _ = unpack_solution(sol)

        u_start = SVector{4,T}(us[1][1:4])
        v_start = SVector{4,T}(us[1][5:8])
        t_start = ts[1]

        u_end = SVector{4,T}(us[end][1:4])
        v_end = SVector{4,T}(us[end][5:8])
        t_end = ts[end]

        GeodesicPoint(sol.prob.p.status, t_start, t_end, u_start, u_end, v_start, v_end)
    end
end

"""
    $(TYPEDSIGNATURES)

Unpacks each point in the solution, similar to [`process_solution`](@ref) but returns an
array of [`GeodesicPoint`](@ref).
"""
function process_solution_full(
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
            GeodesicPoint(sol.prob.p.status, t_start, ti, u_start, ui, v_start, vi)
        end
    end
end

# TODO: GeodesicPath structure for the full geodesic path
# do we want to support this?

function unpack_solution(sol::SciMLBase.AbstractODESolution)
    u = sol.u
    p = sol.prob.p
    t = sol.t
    (u, t, p)
end

unpack_solution(simsol::SciMLBase.AbstractEnsembleSolution) = map(unpack_solution, simsol.u)
