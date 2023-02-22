using Gradus
using Test

# chosen to be close to naked singularity
m = JohannsenPsaltisMetric(M = 1.0, a = 0.8831, ϵ3 = 0.4)
u = SVector(0.0, 1000.0, π / 2, 0.0)

_, _, img = rendergeodesics(m, u, 2000.0, image_width = 100, image_height = 100, fov = 6.0)
# default chart
@test sum(filter(!isnan, img)) ≈ 2.7725698197810156e6 rtol = 0.0001

_, _, img = rendergeodesics(
    m,
    u,
    2000.0,
    image_width = 100,
    image_height = 100,
    fov = 6.0,
    chart = event_horizon_chart(m, closest_approach = 1.001, resolution = 100),
)
# chart that maps event horizon shape better
@test sum(filter(!isnan, img)) ≈ 2.762717446130311e6 rtol = 0.0001
