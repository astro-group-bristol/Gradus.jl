using Test
using Gradus
using StaticArrays

include("../utils.jl")

m = KerrMetric(M = 1.0, a = 0.6)
u = @SVector [0.0, 1000.0, deg2rad(60), 0.0]
d = GeometricThinDisc(0.0, 250.0, π / 2)

x, y

rtype = returntype(lineprofile, m, u, d)
@test isconcretetype(rtype) || "$(rtype)"
