using Gradus
using BenchmarkTools

tracing_suite = BenchmarkGroup(["tracing"])

m = KerrMetric(M = 1.0, a = 0.0)
u = SVector(0.0, 1000.0, deg2rad(75.0), 0.0)
v = map_impact_parameters(m, u, 4.0, 0.0)
d = GeometricThinDisc(0.0, 1000.0, π / 2)
time_span = (0.0, 2000.0)

# single position, single velocity

tracing_suite["single-geodesic"] = BenchmarkGroup(["tracing"])
tracing_suite["single-geodesic"]["no-disc"] =
    @benchmarkable tracegeodesics($m, $u, $v, $time_span)
tracing_suite["single-geodesic"]["no-disc-no-save"] =
    @benchmarkable tracegeodesics($m, $u, $v, $time_span; save_on = false)
tracing_suite["single-geodesic"]["with-disc"] =
    @benchmarkable tracegeodesics($m, $u, $v, $d, $time_span)
tracing_suite["single-geodesic"]["with-disc-no-save"] =
    @benchmarkable tracegeodesics($m, $u, $v, $d, $time_span; save_on = false)

# many positions, many velocities

vs = vec([
    map_impact_parameters(m, u, a, b) for a in range(-11.0, 11.0, 128) for
    b in range(-11.0, 11.0, 128)
])
us = fill(u, length(vs))

tracing_suite["many-geodesic"] = BenchmarkGroup(["tracing"])
tracing_suite["many-geodesic"]["no-disc"] =
    @benchmarkable tracegeodesics($m, $us, $vs, $time_span)
tracing_suite["many-geodesic"]["no-disc-no-save"] =
    @benchmarkable tracegeodesics($m, $us, $vs, $time_span; save_on = false)
tracing_suite["many-geodesic"]["with-disc"] =
    @benchmarkable tracegeodesics($m, $us, $vs, $d, $time_span)
tracing_suite["many-geodesic"]["with-disc-no-save"] =
    @benchmarkable tracegeodesics($m, $us, $vs, $d, $time_span; save_on = false)
tracing_suite["many-geodesic"]["with-disc-no-save-ensemble-endpoint"] =
    @benchmarkable tracegeodesics(
        $m,
        $us,
        $vs,
        $d,
        $time_span;
        save_on = false,
        ensemble = Gradus.EnsembleEndpointThreads(),
    )
tracing_suite["many-geodesic"]["no-disc-no-save-ensemble-endpoint"] =
    @benchmarkable tracegeodesics(
        $m,
        $us,
        $vs,
        $time_span;
        save_on = false,
        ensemble = Gradus.EnsembleEndpointThreads(),
    )
