module FirstOrderMethods

using Accessors
using Parameters
using DocStringExtensions
using StaticArrays

using SciMLBase
using DiffEqCallbacks

# for doc bindings
import ..Gradus

import ..GradusBase:
    AbstractMetricParams,
    inner_radius,
    AbstractGeodesicPoint,
    unpack_solution,
    SciMLBase,
    getgeodesicpoint

import ..GeodesicTracer:
    DiscreteCallback,
    terminate!,
    integrator_problem,
    metric_callback,
    create_callback_set,
    constrain,
    alpha_beta_to_vel,
    ensure_chart_callback

"""
    AbstractFirstOrderMetricParams{T} <: AbstractMetricParams{T}

Abstract type for metrics using the 1st-order integration method. The 1st-order methods
reuse the velocity vector as a parameter vector, where only element `vel[2]` and `vel[3]`
are used, and are local observer ratios ``\\sin \\Theta`` and ``\\sin \\Phi`` respectively.

Require implementation of
- [`inner_radius`](@ref)
- [`constrain`](@ref)
- [`FirstOrderMethods.four_velocity`](@ref)
- [`FirstOrderMethods.calc_lq`](@ref)
- [`FirstOrderMethods.Vr`](@ref)
- [`FirstOrderMethods.Vθ`](@ref)
- [`alpha_beta_to_vel`](@ref)
"""
abstract type AbstractFirstOrderMetricParams{T} <: AbstractMetricParams{T} end

include("implementation.jl")
include("callbacks.jl")

export AbstractFirstOrderMetricParams, FirstOrderGeodesicPoint, BoyerLindquistFO

end # module
