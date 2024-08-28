using Test
using Gradus

m = KerrMetric(1.0, 0.998)
lp = LampPostModel(h = 10.0)
d = ThinDisc(0.0, 1000.0)

sols = @time tracegeodesics(
    m,
    lp,
    d,
    (0.0, 2000.0),
    n_samples = 10_000,
    sampler = EvenSampler(BothHemispheres(), RandomGenerator()),
    save_on = false,
)
gps = filter(i -> status(i) == StatusCodes.IntersectedWithGeometry, unpack_solution(sols))

# ensure error is thrown if not sorted
@test_throws "" RadialDiscProfile(m, lp, gps; N = 100, grid = Gradus.Grids.GeometricGrid())

sort!(gps; by = i -> i.x[2])
prof = RadialDiscProfile(m, lp, gps; N = 100, grid = Gradus.Grids.GeometricGrid())

bins, emiss = get_emissivity(prof)

# get the points and just compute a generic check
@test sum(emiss) ≈ 2025.04 atol = 1e-2
