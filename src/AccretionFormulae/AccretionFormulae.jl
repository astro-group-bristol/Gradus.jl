module AccretionFormulae

import ..Gradus
import ..GradusBase: AbstractMetricParams, metric
import ..AccretionGeometry: CircularOrbits
using ..FirstOrderMethods: FirstOrderGeodesicPoint
using ..Rendering: PointFunction

using Optim: optimize, minimizer, GoldenSection, Brent
using DocStringExtensions
using StaticArrays
using Interpolations
using Tullio: @tullio

include("redshift.jl")
include("orbit-discovery.jl")
include("orbit-interpolations.jl")

export solve_equitorial_circular_orbit,
    trace_equitorial_circular_orbit,
    CircularOrbits,
    interpolate_plunging_velocities,
    PlungingInterpolation,
    interpolate_redshift

end # module
