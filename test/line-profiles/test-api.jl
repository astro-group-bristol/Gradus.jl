using Test
using Gradus
using StaticArrays

include("../utils.jl")

m = KerrMetric(M = 1.0, a = 0.6)
u = @SVector [0.0, 1000.0, deg2rad(60), 0.0]
d = GeometricThinDisc(0.0, 250.0, π / 2)

rtype = returntype(lineprofile, m, u, d)
@test isconcretetype(rtype) || "$(rtype)"

prof = emissivity_profile(m, d, LampPostModel(), n_samples = 100)

x, y = lineprofile(m, u, d, prof; numrₑ = 3, N = 20, algorithm = CunninghamLineProfile())
@test sum(y) ≈ 1 atol = 1e-4

x, y = lineprofile(
    m,
    u,
    d,
    prof;
    plane = PolarPlane(GeometricGrid(); Nr = 10, Nθ = 10, r_max = 90.0),
    algorithm = BinnedLineProfile(),
)
@test sum(y) ≈ 1 atol = 1e-4
