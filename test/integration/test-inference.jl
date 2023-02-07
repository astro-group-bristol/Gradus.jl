using Test
using Gradus
using StaticArrays

include("../utils.jl")

m = KerrMetric()
u = @SVector [0.0, 1e3, π / 2, 0.0]
v = @SVector [0.0, 1.0, 0.0, 0.0]

# single vector for u and v
rtype = returntype(tracegeodesics, m, u, v, (0.0, 2000.0))
@test isconcretetype(rtype)

# multiple vectors for u and v
u = [u, u]
v = [v, v]
rtype = returntype(tracegeodesics, m, u, v, (0.0, 2000.0))
@test isconcretetype(rtype)

# v is a velocity function
u = @SVector [0.0, 1e3, π / 2, 0.0]
vfunc = (i) -> @SVector [0.0, 1.0, 0.0, 0.0]
rtype = returntype((args...) -> tracegeodesics(args...; trajectories = 10), m, u, vfunc, (0.0, 2000.0))
@test isconcretetype(rtype)

# v is a plane
plane = PolarPlane(LinearGrid())
rtype = returntype(tracegeodesics, m, u, plane, (0.0, 2000.0))
@test isconcretetype(rtype)
