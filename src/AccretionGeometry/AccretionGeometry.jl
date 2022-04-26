module AccretionGeometry

import ..GeodesicTracer: tracegeodesics, DiscreteCallback, terminate!
import ..Rendering: rendergeodesics, prerendergeodesics
import Base: in
import RecursiveArrayTools: ArrayPartition
import GeometryBasics
import LinearAlgebra: ×, ⋅

using ..GradusBase
using StaticArrays

include("geometry.jl")
include("meshes.jl")
include("intersections.jl")
include("discs.jl")

function tracegeodesics(
    m::AbstractMetricParams{T},
    init_positions,
    init_velocities,
    accretion_geometry,
    time_domain::Tuple{T,T};
    callback = nothing,
    kwargs...,
) where {T}
    cbs = add_collision_callback(callback, accretion_geometry)
    tracegeodesics(
        m,
        init_positions,
        init_velocities,
        time_domain;
        callback = cbs,
        kwargs...,
    )
end

function rendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    accretion_geometry,
    max_time::T;
    callback = nothing,
    kwargs...,
) where {T}
    cbs = add_collision_callback(callback, accretion_geometry)
    rendergeodesics(m, init_pos, max_time; callback = cbs, kwargs...)
end

function prerendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    accretion_geometry,
    max_time::T;
    callback = nothing,
    kwargs...,
) where {T}
    cbs = add_collision_callback(callback, accretion_geometry)
    prerendergeodesics(m, init_pos, max_time; callback = cbs, kwargs...)
end

add_collision_callback(::Nothing, accretion_geometry) =
    build_collision_callback(accretion_geometry)
add_collision_callback(callback::Base.AbstractVecOrTuple, accretion_geometry) =
    (callback..., build_collision_callback(accretion_geometry))
add_collision_callback(callback, accretion_geometry) =
    (callback, build_collision_callback(accretion_geometry))

function build_collision_callback(geometry::AbstractAccretionGeometry{T}) where {T}
    DiscreteCallback(collision_callback(geometry), i -> terminate!(i, :Intersected))
end

collision_callback(m::AbstractAccretionGeometry{T}) where {T} =
    (u, λ, integrator) -> intersects_geometry(m, line_element(u, integrator))


export AbstractAccretionGeometry,
    AbstractAccretionDisc, MeshAccretionGeometry, GeometricThinDiscm, in_polygon, get_area

end # module
