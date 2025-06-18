using Test
using Gradus

m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 10_000, deg2rad(75), 0.0)
d = ShakuraSunyaev(m)

tf = cunningham_transfer_function(m, x, d, 3.0; β₀ = 1.0)

total = sum(filter(!isnan, tf.f))
@test total ≈ 11.729897381439205 atol = 1e-4

m = KerrMetric(1.0, 0.2)
x = SVector(0.0, 10_000, deg2rad(20), 0.0)
d = ShakuraSunyaev(m; eddington_ratio = 0.2)

tf = cunningham_transfer_function(m, x, d, 5.469668466100368; β₀ = 1.0)
total = sum(filter(!isnan, tf.f))
@test total ≈ 19.321053039396688 atol = 1e-2

# the transfer function here is pretty horrible as it's almost impossible to actually 
# see; this is a test to make sure it doesn't error
# an offset to the isco of 4-e2 resolves this, but that's quite a lot
tf = cunningham_transfer_function(m, x, d, Gradus.isco(m) + 1e-2; β₀ = 1.0)
