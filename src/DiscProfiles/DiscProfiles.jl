module DiscProfiles

using Parameters

import StaticArrays: @SVector, SVector


import ..GradusBase:
    AbstractMetricParams, GeodesicPoint, vector_to_local_sky, getendpoint, inner_radius
import ..GeodesicTracer: tracegeodesics
import ..AccretionGeometry:
    AbstractAccretionGeometry,
    AbstractAccretionDisc,
    to_cartesian,
    inpolygon,
    getarea,
    getcycliclines,
    getpoints

import SciMLBase
import VoronoiCells
import ThreadsX
import GeometryBasics

import Base

include("sky-geometry.jl")
include("corona-models.jl")
include("disc-profiles.jl")
include("transfer-functions.jl")

function tracegeodesics(
    m::AbstractMetricParams{T},
    model::AbstractCoronaModel{T},
    time_domain::Tuple{T,T};
    n_samples = 1024,
    sampler = WeierstrassSampler(res = 100.0),
    kwargs...,
) where {T}
    us = sample_position(m, model, n_samples)
    vs = sample_velocity(m, model, sampler, us, n_samples)
    tracegeodesics(m, us, vs, time_domain; kwargs...)
end


function tracegeodesics(
    m::AbstractMetricParams{T},
    model::AbstractCoronaModel{T},
    d::AbstractAccretionGeometry{T},
    time_domain::Tuple{T,T},
    ;
    n_samples = 1024,
    sampler = WeierstrassSampler(res = 100.0),
    kwargs...,
) where {T}
    us = sample_position(m, model, n_samples)
    vs = sample_velocity(m, model, sampler, us, n_samples)
    tracegeodesics(m, us, vs, d, time_domain; kwargs...)
end


function renderprofile(
    m::AbstractMetricParams{T},
    model::AbstractCoronaModel{T},
    d::AbstractAccretionGeometry{T},
    max_time::T;
    n_samples = 1024,
    sampler = WeierstrassSampler(res = 100.0),
    kwargs...,
) where {T}
    __renderprofile(m, model, d, n_samples, (0.0, max_time); kwargs...)
end

export AbstractCoronaModel,
    LampPostModel,
    renderprofile,
    LowerHemisphere,
    BothHemispheres,
    EvenSampler,
    WeierstrassSampler,
    RandomGenerator,
    GoldenSpiralGenerator,
    VoronoiDiscProfile,
    findindex,
    getareas,
    bin_transfer_function


end # module
