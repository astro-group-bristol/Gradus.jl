using Test
using Gradus
using StaticArrays

m = BoyerLindquistAD()
u = @SVector [0.0, 1e3, π / 2, 0.0]
v = @SVector [0.0, 1.0, 0.0, 0.0]

# single vector for u and v
info = code_typed(tracegeodesics, (typeof(m), typeof(u), typeof(v), Tuple{Float64,Float64}))
info = first(info)
@test isconcretetype(info.second)

# multiple vectors for u and v
u = [u, u]
v = [v, v]
info = code_typed(tracegeodesics, (typeof(m), typeof(u), typeof(v), Tuple{Float64,Float64}))
info = first(info)
@test isconcretetype(info.second)

# v is a velocity function
u = @SVector [0.0, 1e3, π / 2, 0.0]
vfunc = (i) -> @SVector [0.0, 1.0, 0.0, 0.0]
info = code_typed(
    (args...) -> tracegeodesics(args...; trajectories = 10),
    (typeof(m), typeof(u), typeof(vfunc), Tuple{Float64,Float64}),
)
info = first(info)
@test isconcretetype(info.second)

# v is a plane
plane = PolarPlane(LinearGrid())
info = code_typed(
    tracegeodesics,
    (typeof(m), typeof(u), typeof(plane), Tuple{Float64,Float64}),
)
info = first(info)
@test isconcretetype(info.second)
