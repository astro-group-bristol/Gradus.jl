using Test, Gradus

r, θ = Gradus.oblate_spheroid_to_spherical(1.02, 1.113, 0.998)

@test r ≈ 1.3872 atol = 1e-3
@test θ ≈ acos(0.8023) atol = 1e-3
