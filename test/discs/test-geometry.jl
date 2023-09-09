using Gradus, Test

m = KerrMetric(1.0, 0.9)
d = ShakuraSunyaev(m)
n = Gradus._cartesian_tangent_vector(d, 2.6)
@test n ≈ [0.689693957000099, 0.0, 0.724100991352412] atol = 1e-5

n = Gradus._cartesian_tangent_vector(d, 6.6)
@test n ≈ [0.9679192396299138, 0.0, 0.2512615083021063] atol = 1e-5

n = Gradus._cartesian_tangent_vector(d, 1.0)
@test n ≈ [1, 0, 0] atol = 1e-5
