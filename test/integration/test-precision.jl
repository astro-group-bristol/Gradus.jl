using Test
using Gradus


m = KerrMetric(M = 1.0, a = 1.0)
u = SVector(0.0, 1000.0, π / 2, 0.0)

# target position
target = SVector(10.0, 0.005, 0.0)

α, β, accuracy = Gradus.impact_parameters_for_target(target, m, u)
@test α ≈ -0.004013630261097743 rtol = 1e-3
@test β ≈ 10.969606493445841 rtol = 1e-3
@test accuracy ≈ 0.004587323209289997 rtol = 1e-3

# new target
target = SVector(10.0, deg2rad(40), -π / 4)

α, β, accuracy = Gradus.impact_parameters_for_target(target, m, u)
@test α ≈ 4.848373364532467 rtol = 1e-3
@test β ≈ 8.02066263349774 rtol = 1e-3
@test accuracy ≈ 0.0017037999175873982 rtol = 1e-3

# test geodesic point interface
α, β, gp, accuracy = Gradus.optimize_for_target(target, m, u)
@test gp.x[1] ≈ 1005.2700874611182 rtol = 1e-3
@test α ≈ 4.848373364532467 rtol = 1e-3
@test β ≈ 8.02066263349774 rtol = 1e-3
@test accuracy ≈ 0.0017037999175873982 rtol = 1e-3
