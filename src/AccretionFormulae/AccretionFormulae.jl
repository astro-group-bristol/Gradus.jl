module AccretionFormulae

import ..Gradus
import ..GradusBase: AbstractMetricParams, metric
using ..FirstOrderMethods: FirstOrderGeodesicPoint
using ..Rendering: PointFunction

using Optim: optimize, minimizer, GoldenSection, Brent
using DocStringExtensions
using StaticArrays
using Interpolations

include("redshift.jl")
include("orbit-discovery.jl")
include("orbit-interpolations.jl")

export solve_equitorial_circular_orbit,
    trace_equitorial_circular_orbit,
    CircularOrbits,
    interpolate_plunging_velocities,
    PlungingInterpolation

end # module
