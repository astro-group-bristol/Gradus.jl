using Test
using Gradus

m = KerrMetric(1.0, 0.2)
d = PolishDoughnut(m, 12.0, 0.21);
d.inner_radius;

r = collect(range(10.0, 15.0, 200))
h = map(x -> Gradus.r_cross_section(d, x), r)

# fingerprint the cross section map
@test sum(h) ≈ 219.97440610254944 atol = 1e-5
