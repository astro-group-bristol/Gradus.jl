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

include("carter-method-bl-impl.jl")
include("carter-method-bl-interface.jl")
include("implementation.jl")
include("callbacks.jl")

"""
$(TYPEDSIGNATURES)

Radius of marginally stable orbit.

From Bardeen et al. (1972) eq. (2.21):

```math
r_\\text{ms} = M \\left\\{ 3 + Z_2 \\pm \\sqrt{(3 - Z_1)(3 + Z_1 + 2 Z_2)} \\right\\}.
```

The choice of ``\\pm`` is chosen by the sign of ``a``.
"""
rms(M, a, ±) = M * (3 + Z₂(M, a) ± √((3 - Z₁(M, a)) * (3 + Z₁(M, a) + 2 * Z₂(M, a))))
rms(M, a) = a > 0.0 ? rms(M, a, -) : rms(M, a, +)

end # module
