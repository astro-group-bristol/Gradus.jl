using Test
using Gradus

m = KerrMetric(1.0, 0.998)
d = GeometricThinDisc(0.0, 100.0, π / 2)

model = LampPostModel(h = 10.0)

# test usual case
em_prof1 = Gradus.emissivity_profile(
    m,
    d,
    model,
)

# try explicitly setting a spectrum
spectrum = Gradus.PowerLawSpectrum(3.0)

# now with spectrum as an argument
em_prof2 = Gradus.emissivity_profile(
    m,
    d,
    model,
    spectrum,
)

