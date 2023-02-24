using Test
using Gradus


m = KerrMetric(M = 1.0, a = 1.0)
u = SVector(0.0, 1000.0, π / 2, 0.0)

# target position
target = SVector(10.0, 0.001, 0.0)

α, β, accuracy = Gradus.impact_parameters_for_target(target, m, u)
@test α ≈ -0.001950232871615813 rtol = 1e-3
@test β ≈ 10.969606493445841 rtol = 1e-3
@test accuracy ≈ 0.00870713834731178 rtol = 1e-3

# new target
target = SVector(10.0, deg2rad(40), -π / 4)

α, β, accuracy = Gradus.impact_parameters_for_target(target, m, u)
@test α ≈ 4.857623042719467 rtol = 1e-3
@test β ≈ 8.025434644810364 rtol = 1e-3
@test accuracy ≈ 0.0021553188931118686 rtol = 1e-3
