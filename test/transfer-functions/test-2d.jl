using Test
using Gradus
using StaticArrays

m = KerrMetric(M = 1.0, a = 0.998)
u = @SVector [0.0, 1e6, deg2rad(30), 0.0]

plane = PolarPlane(GeometricGrid(), Nr = 20, Nθ = 20)

d = GeometricThinDisc(Gradus.isco(m), 500.0, π / 2)
model = LampPostModel(h = 10.0, θ = deg2rad(0.0001))

# calculate transfer functions
tf = lagtransfer(
    m,
    u,
    d,
    model;
    plane = plane,
    callback = domain_upper_hemisphere(),
    n_samples = 100,
    sampler = EvenSampler(domain = BothHemispheres(), generator = GoldenSpiralGenerator()),
)

# check number of intersection points
@test length(tf.observer_to_disc) == 311
@test length(tf.emissivity_profile.geodesic_points) == 58

# ensure binning works as expected
t, E, f = binflux(tf, N_t = 100, N_E = 100)

fluxsum = sum(filter(!isnan, f))
@test fluxsum ≈ 8.5152 atol = 1e-2
