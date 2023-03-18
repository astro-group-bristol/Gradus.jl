using Test
using Gradus

# test that tracing with charged test-particles actually
# produces different results (/ is being executed correctly)
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
        fov = 2.5,
        q = q,
        kwargs...,
    )
    fingerprint = sum(filter(!isnan, img))
    fingerprint
end

@test test_tracer(m, u, 0.0) ≈ 450413.5565857728 rtol = 1e-3
@test test_tracer(m, u, 1.0) ≈ 263050.89536876464 rtol = 1e-3
@test test_tracer(m, u, -1.0) ≈ 652669.0871285137 rtol = 1e-3

ensemble = Gradus.EnsembleEndpointThreads()
@test test_tracer(m, u, 0.0; ensemble = ensemble) ≈ 450413.5565857728 rtol = 1e-3
@test test_tracer(m, u, 1.0; ensemble = ensemble) ≈ 263050.89536876464 rtol = 1e-3
@test test_tracer(m, u, -1.0; ensemble = ensemble) ≈ 652669.0871285137 rtol = 1e-3
