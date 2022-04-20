module DiscProfiles

using Parameters

import StaticArrays: @SVector


import ..GradusBase: AbstractMetricParams, vector_to_local_sky
import ..GeodesicTracer: tracegeodesics
import ..AccretionGeometry: AbstractAccretionGeometry

import GeometricalPredicates: Point2D, getx, gety, geta, getb, getc

include("sky-geometry.jl")
include("corona-models.jl")

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
    time_domain::Tuple{T,T},
    d::AbstractAccretionGeometry{T},
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
    max_time::T,
    d::AbstractAccretionGeometry{T};
    n_samples = 2048,
    sampler = WeierstrassSampler(res = 100.0),
    kwargs...,
) where {T}
    __renderprofile(m, model, d, n_samples, (0.0, max_time); kwargs...)
end

export AbstractCoronaModel,
    LampPostModel,
    tracegeodesics,
    renderprofile,
    LowerHemisphere,
    BothHemispheres,
    EvenSampler,
    WeierstrassSampler,
    RandomGenerator,
    GoldenSpiralGenerator


end # module
