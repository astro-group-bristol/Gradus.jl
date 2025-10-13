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
     1.4346387869787864,
     3.0822515234888774,
     1.7923604648828981,
     0.6016959946033558,
     0.11910008907351012,
     0.017392602799041507,
     0.0023309504405384547,
     0.0003139154565507922,
     3.665392374360994e-5,
     1.2069687133228597e-6,
] rtol = 1e-2
