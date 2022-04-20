module FirstOrderMethods

using Accessors
using Parameters
using DocStringExtensions
using StaticArrays

using SciMLBase
using DiffEqCallbacks

import ..GradusBase:
    AbstractMetricParams,
    inner_radius,
    AbstractGeodesicPoint,
    get_endpoint,
    geodesic_point_type,
    unpack_solution,
    SciMLBase

import ..GeodesicTracer:
    DiscreteCallback,
    terminate!,
    integrator_problem,
    metric_callback,
    create_callback_set,
    constrain,
    alpha_beta_to_vel

abstract type AbstractFirstOrderMetricParams{T} <: AbstractMetricParams{T} end

include("implementation.jl")
include("callbacks.jl")

export AbstractFirstOrderMetricParams, FirstOrderGeodesicPoint, BoyerLindquistFO

end # module
