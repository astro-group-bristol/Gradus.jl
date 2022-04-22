module AccretionFormulae

import ..Gradus
import ..GradusBase: AbstractMetricParams, metric

using ..FirstOrderMethods: FirstOrderGeodesicPoint
using ..Rendering: PointFunction

using Optim: optimize, minimizer, GoldenSection, Brent
using DocStringExtensions
using StaticArrays

include("redshift.jl")
include("orbit-discovery.jl")

const redshift = PointFunction(_redshift_guard)

export solve_equitorial_circular_orbit, trace_equitorial_circular_orbit, redshift

end # module
