using Test
using Gradus

m = KerrMetric(1.0, 0.998)
model = LampPostModel(h = 10.0)
d = ThinDisc(0.0, 500.0)

# using point source angular method
profile = emissivity_profile(m, d, model; n_samples = 20)

r, em = profile.radii, profile.ε
@test em ≈ [
    0.0029464479567890534,
    0.0014052519492578114,
    0.0008963679521766861,
    0.0005749351642563003,
    0.0003386885861792927,
    0.0001703542742784169,
    6.482839568020104e-5,
    1.3029008103481133e-5,
    3.432060732289487e-6,
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
    1.4803144540690316,
    3.2088064526759394,
    1.9690622763490204,
    0.8610923964275128,
    0.2343074914871196,
    0.04489102975754837,
    0.007577479383680966,
    0.0012209856249765652,
    0.0002103466340875971,
    3.3564404237996266e-5,
] rtol = 1e-2
