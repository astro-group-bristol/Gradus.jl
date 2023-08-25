using Gradus
using Test

m = KerrMetric(1.0, 0.998)
d = GeometricThinDisc(0.0, 100.0, π / 2)

model0 = LampPostModel(h = 10.0)
model1 = Gradus.BeamedPointSource(10.0, 0.0)

em_prof0 = Gradus.emissivity_profile(m, d, model0, n_samples = 100)

em_prof1 = Gradus.emissivity_profile(m, d, model1, n_samples = 100)

profile0 = Gradus.RadialDiscProfile(em_prof0)

profile1 = Gradus.RadialDiscProfile(em_prof1)

test_radii = range(2, 100, 10) |> collect
@test profile0.f.ε.(test_radii) ≈ profile1.f.ε.(test_radii) rtol = 1e-1
