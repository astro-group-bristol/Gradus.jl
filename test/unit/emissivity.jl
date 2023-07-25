using Test
using Gradus

m = KerrMetric(1.0, 0.998)
model = LampPostModel(h = 10.0)
d = GeometricThinDisc(0.0, 500.0, π / 2)

# using point source angular method
profile = emissivity_profile(m, d, model; n_samples = 20)

r, em = get_emissivity(profile)
@test em ≈ [
    0.03530163295788805,
    0.01618161075710721,
    0.010034247429248367,
    0.006231511729195148,
    0.0035296670918803694,
    0.0016859014503211888,
    0.0005946901244498561,
    0.0001052324320767813,
] atol = 1e-5

# using generic monte-carlo sampling
profile = emissivity_profile(
    m,
    d,
    model;
    n_samples = 1000,
    sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
    N = 10,
)

r, em = get_emissivity(profile)

@test em ≈ [
    1.4536471470344825,
    3.2159234810191997,
    1.9624518676973357,
    0.8654137966466704,
    0.23973494196571282,
    0.046232848740113686,
    0.00789393186371828,
    0.0012974087513837867,
    0.0001872022563608107,
    3.825853506250807e-5,
] atol = 1e-5
