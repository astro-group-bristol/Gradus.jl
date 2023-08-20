using Test
using Gradus

m = KerrMetric(1.0, 0.998)
model = LampPostModel(h = 10.0)
d = GeometricThinDisc(0.0, 500.0, π / 2)

# using point source angular method
profile = emissivity_profile(m, d, model; n_samples = 20)

r, em = get_emissivity(profile)
@test em ≈ [
    0.035608673876360804,
    0.016220783046579302,
    0.01004833455994888,
    0.006238633814155616,
    0.0035332505361847172,
    0.0016872870193954522,
    0.00059491166015569,
    0.00010511998411323895,
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
