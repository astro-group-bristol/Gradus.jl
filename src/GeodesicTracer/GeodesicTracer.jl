module GeodesicTracer

using SciMLBase
using OrdinaryDiffEq
using DiffEqCallbacks

using StaticArrays
using DocStringExtensions
using Parameters

import ..GradusBase: AbstractMetricParams, geodesic_eq, constrain, inner_radius, metric

import ForwardDiff

include("callbacks.jl")
include("problem.jl")
include("tracer.jl")
include("constraints.jl")
include("utility.jl")
include("auto-diff.jl")

"""
    tracegeodesics(
        m::AbstractMetricParams{T},
        position, velocity,
        time_domain::Tuple{T,T}
        ;
        μ = 0.0f0,
        callbacks=Nothing,
        solver=Tsit5(),
        solver_opts...
    )

Integrate a geodesic for metric parameterised by `m`, for some initial positions and velocities.
The positions and velocities may be

  - a single position and velocity in the form of a vector of numbers
  - a collection of positions and velocities, as either a vector of vectors, or as a matrix

The matrix specification reads each corresponding column as the initial position and velocity. When a collection of
positions and velocities is supplied, this method dispatched `EnsembleProblem`, offering `ensemble` as a `solver_opts`,
specifying the ensemble method to use.

`solver_opts` are the common solver options in DifferentialEquations.
"""
function tracegeodesics(
    m::AbstractMetricParams{T},
    position,
    velocity,
    time_domain::Tuple{T,T};
    solver = Tsit5(),
    μ = 0.0,
    closest_approach = 1.01,
    effective_infinity = 1200.0,
    callback = nothing,
    solver_opts...,
) where {T}

    _velocity = if (velocity isa Function) && (eltype(position) === T)
        wrap_constraint(m, position, velocity, μ)
    else
        if eltype(position) !== eltype(velocity)
            error("Position and velocity must have the same element type.")
        end
        constrain_all(m, position, velocity, μ)
    end

    cbs = create_callback_set(m, callback, closest_approach, effective_infinity)

    __tracegeodesics(
        m,
        position,
        _velocity,
        time_domain,
        solver;
        callback = cbs,
        abstol = 1e-9,
        reltol = 1e-9,
        solver_opts...,
    )
end

export tracegeodesics,
    map_impact_parameters,
    AbstractAutoDiffMetricParams,
    AbstractAutoDiffStaticAxisSymmetricParams,
    metric_components,
    metric_jacobian

end # module
