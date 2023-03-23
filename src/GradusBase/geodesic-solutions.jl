abstract type AbstractIntegrationParameters end

update_integration_parameters!(
    p::AbstractIntegrationParameters,
    ::AbstractIntegrationParameters,
) = error("Not implemented for $(typeof(p))")

mutable struct IntegrationParameters <: AbstractIntegrationParameters
    status::StatusCodes.T
end

function update_integration_parameters!(
    p::IntegrationParameters,
    new::IntegrationParameters,
)
    p.status = new.status
    p
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
struct GeodesicPoint{T,A} <: AbstractGeodesicPoint{T}
    "Return code of the integrator for this geodesic."
    status::StatusCodes.T
    "Start affine time"
    λ_min::T
    "End affine time"
    λ_max::T
    "Start four position"
    x_init::SVector{4,T}
    "End four position"
    x::SVector{4,T}
    "Start four velocity"
    v_init::SVector{4,T}
    "End four velocity"
    v::SVector{4,T}
    "Auxillary values (polarisation, intensities, etc.)"
    aux::A
end

function Base.show(io::IO, gp::GeodesicPoint)
    print(io, "GeodesicPoint($(gp.x...))")
end

function Base.show(io::IO, ::MIME"text/plain", gp::GeodesicPoint)
    # text = """GeodesicPoint:
    #   . status : $(gp.status)
    #   . λ_min  : $(gp.λ_min)
    #   . λ_max  : $(gp.λ_max)
    #   . x_init : $(gp.x_init)
    #   . v_init : $(gp.v_init)
    #   . x      : $(gp.x)
    #   . v      : $(gp.v)
    # """
    text =
        "GeodesicPoint:\n" * join(
            (
                rpad("  . $f", 12) * ": $(getproperty(gp, f))" for
                f in fieldnames(typeof(gp))
            ),
            "\n",
        )
    print(io, text)
end

"""
    process_solution([m], sol)

Unpack a solution (`SciMLBase.AbstractODESolution`) as a [`GeodesicPoint`](@ref), optionally specifying
the metric under which quantities are transformed. 

If the solution stores any additional parameters (e.g. intensity in radiative transfer), these will be packed 
into the `aux` field of [`GeodesicPoint`](@ref).

## Example use

```julia
sol = tracegeodesics(m, x, v)
point = process_solution(sol)
```
"""
function process_solution(_, sol::SciMLBase.AbstractODESolution{T}) where {T}
    @inbounds @views begin
        us, ts, _ = unpack_solution(sol)

        u_start = SVector{4,T}(us[1][1:4])
        v_start = SVector{4,T}(us[1][5:8])
        t_start = ts[1]

        u_end = SVector{4,T}(us[end][1:4])
        v_end = SVector{4,T}(us[end][5:8])
        t_end = ts[end]

        GeodesicPoint(
            sol.prob.p.status,
            t_start,
            t_end,
            u_start,
            u_end,
            v_start,
            v_end,
            eltype(us) <: SVector{9} ? us[end][9] : nothing,
        )
    end
end

"""
    $(TYPEDSIGNATURES)

Unpacks each point in the solution, similar to [`process_solution`](@ref) but returns an
array of [`GeodesicPoint`](@ref).
"""
function process_solution_full(
    _::AbstractMetric{T},
    sol::SciMLBase.AbstractODESolution{T},
) where {T}
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

process_solution(gp::AbstractGeodesicPoint) = gp
process_solution(sol::SciMLBase.AbstractODESolution) = process_solution(sol.prob.f.f.m, sol)

# TODO: GeodesicPath structure for the full geodesic path
# do we want to support this?

function unpack_solution(sol::SciMLBase.AbstractODESolution)
    u = sol.u
    p = sol.prob.p
    t = sol.t
    (u, t, p)
end

unpack_solution(simsol::SciMLBase.AbstractEnsembleSolution) = map(unpack_solution, simsol.u)
