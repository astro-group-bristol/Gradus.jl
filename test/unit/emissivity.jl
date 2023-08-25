using Test
using Gradus

m = KerrMetric(1.0, 0.998)
model = LampPostModel(h = 10.0)
d = GeometricThinDisc(0.0, 500.0, π / 2)

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
    1.4795830708196946,
    3.216264925932792,
    1.9626618445135773,
    0.8655046521104081,
    0.23975966603033919,
    0.04623755540947905,
    0.007894729167415742,
    0.0012975391924250703,
    0.0001872210347441914,
    3.8919060688680126e-5,
] rtol = 1e-2
