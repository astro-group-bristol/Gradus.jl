using Test
using Gradus

m = KerrMetric(1.0, 0.998)
model = LampPostModel(h = 10.0)
d = ThinDisc(0.0, 500.0)

# using point source angular method
profile = emissivity_profile(m, d, model; n_samples = 20)

r, em = profile.radii, profile.ε
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

r, em = profile.radii, profile.ε
@test em ≈ [
    1.4796037159131263,
    3.2366973461371877,
    1.9268465534519976,
    0.8551895569098363,
    0.22571555386177453,
    0.04291156396432442,
    0.007185350120798265,
    0.0010608325559791855,
    0.00017033140353248102,
    3.1119048300406384e-5,
] rtol = 1e-2
