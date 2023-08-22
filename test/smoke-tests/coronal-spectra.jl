using Gradus
using Plots
using Test

m = KerrMetric(1.0, 0.998)
d = GeometricThinDisc(0.0, 100.0, π/2)

model = LampPostModel(h = 10.0)

#calculating emissivity profile without defining the spectrum to see if it uses the default Γ value (Γ = 2)
em_prof1 = Gradus.emissivity_profile(
m, 
d, 
model, 
n_samples = 100_000, 
sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
)

#first profile, default photon index Γ
profile1 = Gradus.RadialDiscProfile(em_prof1)

#defining the spectrum
spectrum = PowerLawSpectrum(2.0)

#secondly, passing spectrum as an argument 
em_prof2 = Gradus.emissivity_profile(
m, 
d, 
model, 
spectrum,
n_samples = 100_000, 
sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
)

#profile number 2 with defined Γ = 2
profile2 = Gradus.RadialDiscProfile(em_prof2)

#testing if they are equal
@test profile1.f.ε.u ≈ profile2.f.ε.u
