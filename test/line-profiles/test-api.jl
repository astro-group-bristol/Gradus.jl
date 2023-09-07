using Test
using Gradus
using StaticArrays

include("../utils.jl")

m = KerrMetric(M = 1.0, a = 0.6)
u = @SVector [0.0, 1000.0, deg2rad(60), 0.0]
d = ThinDisc(0.0, 250.0)

rtype = returntype(lineprofile, m, u, d)
@test isconcretetype(rtype) || "$(rtype)"

prof = emissivity_profile(m, d, LampPostModel(), n_samples = 100)

x, y = lineprofile(m, u, d, prof; numrₑ = 3, N = 20, method = TransferFunctionMethod())
@test sum(y) ≈ 1 atol = 1e-4

x, y = lineprofile(
    m,
    u,
    d,
    prof;
    plane = PolarPlane(GeometricGrid(); Nr = 10, Nθ = 10, r_max = 90.0),
    method = BinningMethod(),
)
@test sum(y) ≈ 1 atol = 1e-4
