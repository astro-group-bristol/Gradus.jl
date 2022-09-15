module AccretionFormulae

import ..Gradus
import ..GradusBase: AbstractMetricParams, metric, getgeodesicpoint
import ..GeodesicTracer: map_impact_parameters, tracegeodesics
import ..AccretionGeometry: CircularOrbits, AbstractAccretionDisc, AbstractAccretionGeometry
using ..FirstOrderMethods: FirstOrderGeodesicPoint
using ..Rendering: PointFunction

using LinearAlgebra: det
using Optim: optimize, minimizer, GoldenSection, Brent
using DocStringExtensions
using StaticArrays
using Interpolations
using Tullio: @tullio
using Parameters

import FiniteDifferences
import ThreadsX
import Roots

include("redshift.jl")
include("orbit-discovery.jl")
include("orbit-interpolations.jl")
include("emission-radii.jl")

export solve_equitorial_circular_orbit,
    trace_equitorial_circular_orbit,
    CircularOrbits,
    interpolate_plunging_velocities,
    PlungingInterpolation,
    interpolate_redshift

end # module
