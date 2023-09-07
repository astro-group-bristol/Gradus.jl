using Gradus
using Test

m = KerrMetric(1.0, 0.998)
d = ThinDisc(0.0, 100.0)

model0 = LampPostModel(h = 10.0)
model1 = Gradus.BeamedPointSource(10.0, 0.0)

em_prof0 = Gradus.emissivity_profile(m, d, model0, n_samples = 100)

em_prof1 = Gradus.emissivity_profile(m, d, model1, n_samples = 100)

profile0 = Gradus.RadialDiscProfile(em_prof0)

profile1 = Gradus.RadialDiscProfile(em_prof1)

test_radii = range(2, 100, 10) |> collect
@test emissivity_at(profile0, test_radii) ≈ emissivity_at(profile1, test_radii) rtol = 1e-1
