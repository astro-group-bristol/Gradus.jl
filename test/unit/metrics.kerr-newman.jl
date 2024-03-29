using Test
using Gradus

# test that tracing with charged test-particles actually
# produces different results (/ is being executed correctly)
m = KerrNewmanMetric(M = 1.0, a = 0.6, Q = 0.6)
u = SVector(0.0, 1000.0, π / 2, 0.0)

function test_tracer(m, u, q; kwargs...)
    α, β, img = rendergeodesics(
        m,
        u,
        # max integration time
        2000.0;
        image_width = 40,
        image_height = 40,
        αlims = (-8, 8),
        βlims = (-8, 8),
        q = q,
        kwargs...,
    )
    fingerprint = sum(filter(!isnan, img))
    fingerprint
end
@test test_tracer(m, u, 0.0) ≈ 428809.9681726607 rtol = 1e-3
@test test_tracer(m, u, 1.0) ≈ 253280.6794972752 rtol = 1e-3
@test test_tracer(m, u, -1.0) ≈ 619335.5670363897 rtol = 1e-3

ensemble = Gradus.EnsembleEndpointThreads()
@test test_tracer(m, u, 0.0; ensemble = ensemble) ≈ 428809.9681726607 rtol = 1e-3
@test test_tracer(m, u, 1.0; ensemble = ensemble) ≈ 253280.6794972752 rtol = 1e-3
@test test_tracer(m, u, -1.0; ensemble = ensemble) ≈ 619335.5670363897 rtol = 1e-3

# ensure the circular velocity calculations work
@test CircularOrbits.fourvelocity(m, 20.0; q = 1.0) ≈
      [1.065341126764724, 0.0, 0.0, 0.007652504287280518] rtol = 1e-5
@test CircularOrbits.fourvelocity(m, 20.0; q = -1.0) ≈
      [1.099562453625687, 0.0, 0.0, 0.015091823134051219] rtol = 1e-5
