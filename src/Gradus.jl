module Gradus

import Base: in
using Base.Threads: @threads
using LinearAlgebra: ×, ⋅, norm, det

using DocStringExtensions
using Parameters

using SciMLBase
using OrdinaryDiffEq
using DiffEqCallbacks
using StaticArrays

using Accessors: @set
using Tullio: @tullio

import ThreadsX
import ForwardDiff
import GeometryBasics

include("GradusBase/GradusBase.jl")
using .GradusBase
import .GradusBase: AbstractGeodesicPoint, GeodesicPoint

include("tracing/tracing.jl")
include("tracing/constraints.jl")
include("tracing/callbacks.jl")
include("tracing/utility.jl")

include("tracing/method-implementations/auto-diff.jl")

include("rendering/cache.jl")
include("rendering/rendering.jl")
include("rendering/utility.jl")

include("tracing/method-implementations/first-order.jl")

include("point-functions/point-functions.jl")

include("orbits/circular-orbits.jl")
include("orbits/orbit-discovery.jl")
include("orbits/orbit-interpolations.jl")

include("accretion-geometry/geometry.jl")
include("accretion-geometry/intersections.jl")
include("accretion-geometry/discs.jl")
include("accretion-geometry/meshes.jl")
include("accretion-geometry/bootstrap.jl")

include("orbits/emission-radii.jl")

include("corona-to-disc/sky-geometry.jl")
include("corona-to-disc/corona-models.jl")
include("corona-to-disc/disc-profiles.jl")
include("corona-to-disc/transfer-functions.jl")

include("metrics/boyer-lindquist-ad.jl")
include("metrics/boyer-lindquist-fo.jl")
include("metrics/johannsen-ad.jl")
include("metrics/johannsen-psaltis-ad.jl")
include("metrics/morris-thorne-ad.jl")
include("metrics/kerr-refractive-ad.jl")
include("metrics/dilaton-axion-ad.jl")

include("special-radii.jl")
include("redshift.jl")

include("point-functions/const-point-functions.jl")

end # module
