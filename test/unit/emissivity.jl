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
    1.4531056469477812,
    3.216028944655333,
    1.9625919375522862,
    0.8654885033740061,
    0.2397572736938699,
    0.04623730854934404,
    0.007894706545002786,
    0.0012975372018286873,
    0.0001872208766790682,
    3.826234997699e-5,
] atol = 1e-5
