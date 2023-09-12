using Test
using Gradus

m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 10_000, deg2rad(75), 0.0)
d = ShakuraSunyaev(m)

tf = cunningham_transfer_function(m, x, d, 3.0; β₀ = 1.0)

total = sum(filter(!isnan, tf.f))
@test total ≈ 12.27895052104429 atol = 1e-4
