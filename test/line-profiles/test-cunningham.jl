using Test
using Gradus
using StaticArrays

m = KerrMetric(M = 1.0, a = 0.6)
u = @SVector [0.0, 1000.0, deg2rad(60), 0.0]
d = GeometricThinDisc(0.0, 250.0, π / 2)

bins = collect(range(0.1, 1.3, 100))
x, y =
    lineprofile(bins, (r) -> r^(-3), m, u, d, CunninghamLineProfile(); N = 40, numrₑ = 30)

# should be around .3 - .4
g_low = x[findfirst(>(0), y)]
@test g_low ≈ 0.355 atol = 0.05

# should be 1.15 - 1.25
g_high = x[end-findfirst(>(0), reverse(y))]
@test g_high ≈ 1.2 atol = 0.05

# area under the curve
@test sum(y) ≈ 1.0

# test for other metrics
m = JohannsenPsaltisMetric(M = 1.0, a = 0.6, ϵ3 = 2.0)
x, y =
    lineprofile(bins, (r) -> r^(-3), m, u, d, CunninghamLineProfile(); N = 40, numrₑ = 30)

# should be around 0.27
g_low = x[findfirst(>(0), y)]
@test g_low ≈ 0.27 atol = 0.05

# should be 1.15 - 1.25
g_high = x[end-findfirst(>(0), reverse(y))]
@test g_high ≈ 1.2 atol = 0.05

# area under the curve
@test sum(y) ≈ 1.0
