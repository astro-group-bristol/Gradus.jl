using Test
using Gradus

m = KerrMetric(M = 1.0, a = 0.998)
x = SVector(0.0, 1e6, deg2rad(30), 0.0)

plane = PolarPlane(GeometricGrid(), Nr = 20, Nθ = 20)

d = ThinDisc(Gradus.isco(m), 500.0)
model = LampPostModel(h = 10.0, θ = deg2rad(0.0001))

# calculate transfer functions
tf = lagtransfer(
    m,
    x,
    d,
    model;
    plane = plane,
    callback = domain_upper_hemisphere(),
    n_samples = 100,
    sampler = EvenSampler(domain = BothHemispheres(), generator = GoldenSpiralGenerator()),
)

# check number of intersection points
@test length(tf.observer_to_disc) == 335
@test length(tf.emissivity_profile.geodesic_points) == 58

# ensure binning works as expected
t, E, f = binflux(tf, N_t = 100, N_E = 100)

fluxsum = sum(filter(!isnan, f))
@test fluxsum ≈ 5.037953486417272 atol = 1e-2

# test semi-analytic method

prof = emissivity_profile(
    m,
    d,
    model;
    n_samples = 5_000,
    sampler = EvenSampler(domain = BothHemispheres(), generator = GoldenSpiralGenerator()),
)

radii = Gradus.Grids._inverse_grid(Gradus.isco(m), 100.0, 5)
d = GeometricThinDisc(0.0, 500.0, π / 2)
itb = @time Gradus.interpolated_transfer_branches(m, x, d, radii)

bins = collect(range(0.0, 1.5, 100))
tbins = collect(range(0, 150.0, 100))

flux = Gradus.integrate_lagtransfer(prof, itb, radii, bins, tbins; t0 = x[2], Nr = 1000)

fluxsum = sum(flux)
# normalisation is missing, so number is big
@test fluxsum ≈ 1.0 atol = 1e-2
