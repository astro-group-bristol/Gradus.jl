
"""
    struct GeodesicPoint <: AbstractGeodesicPoint

Fields:
- `status::StatusCodes.T`: Return code of the integrator for this geodesic.
- `λ_min::T`: Start affine time
- `λ_max::T`: End affine time
- `x_init::SVector{4,T}`: Start four position
- `x::SVector{4,T}`: End four position
- `v_init::SVector{4,T}`: Start four velocity
- `v::SVector{4,T}`: End four velocity
- `aux::A`: Auxillary values (polarisation, intensities, etc.)
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
    unpack_solution([m], sol)

Unpack a solution (`SciMLBase.AbstractODESolution`) as a [`GeodesicPoint`](@ref), optionally specifying
the metric under which quantities are transformed. 

If the solution stores any additional parameters (e.g. intensity in radiative transfer), these will be packed 
into the `aux` field of [`GeodesicPoint`](@ref).

## Example use

```julia
sol = tracegeodesics(m, x, v)
point = unpack_solution(sol)
```
"""
function unpack_solution(::AbstractMetric, sol::SciMLBase.AbstractODESolution{T}) where {T}
    @inbounds @views begin
        us, ts, _ = extract_fields(sol)

        x_init = SVector{4,T}(us[1][1:4])
        v_init = SVector{4,T}(us[1][5:8])
        t_init = ts[1]

        x = SVector{4,T}(us[end][1:4])
        v = SVector{4,T}(us[end][5:8])
        t = ts[end]

        # get the auxiliary values if we have any
        aux = unpack_auxiliary(us[end])

        GeodesicPoint(get_status_code(sol.prob.p), t_init, t, x_init, x, v_init, v, aux)
    end
end

unpack_auxiliary(::SVector{8}) = nothing
function unpack_auxiliary(u::SVector{N,T}) where {N,T}
    @assert N > 8
    SVector{N - 8,T}(u[9:end])
end

function extract_fields(sol)
    u = sol.u
    p = sol.prob.p
    t = sol.t
    (u, t, p)
end

"""
    unpack_solution_full

Unpacks each point in the solution, similar to [`unpack_solution`](@ref) but returns an
array of [`GeodesicPoint`](@ref).
"""
function unpack_solution_full(
    _::AbstractMetric{T},
    sol::SciMLBase.AbstractODESolution{T},
) where {T}
    us, ts, _ = extract_fields(sol)
    @inbounds @views begin
        u_start = SVector{4,T}(us[1][1:4])
        v_start = SVector{4,T}(us[1][5:8])
        t_start = ts[1]
        map(eachindex(us)) do i
            ui = SVector{4,T}(us[i][1:4])
            vi = SVector{4,T}(us[i][5:8])
            ti = ts[i]
            GeodesicPoint(
                get_status_code(sol.prob.p),
                t_start,
                ti,
                u_start,
                ui,
                v_start,
                vi,
            )
        end
    end
end

unpack_solution(gp::AbstractGeodesicPoint) = gp
unpack_solution(sol::SciMLBase.AbstractODESolution) = unpack_solution(sol.prob.f.f.m, sol)
unpack_solution(simsol::SciMLBase.AbstractEnsembleSolution) = map(unpack_solution, simsol.u)
