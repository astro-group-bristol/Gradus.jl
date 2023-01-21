using Test
using Gradus
using StaticArrays

# test known radii
m = KerrMetric(M = 1.0, a = 0.0)
@test Gradus.isco(m) == 6.0

m = KerrMetric(M = 1.0, a = 1.0)
@test Gradus.isco(m) == 1.0

# test event horizon code
m = KerrMetric(M = 1.0, a = 0.0)
rs, _ = event_horizon(m)
@test rs ≈ fill(2.0, length(rs))

m = KerrMetric(M = 1.0, a = 1.0)
rs, _ = event_horizon(m)
@test rs ≈ fill(1.0, length(rs))

# make sure works for other metrics too
for M in [JohannsenMetric, JohannsenPsaltisMetric]
    @test Gradus.isco(M()) > 1.0
    _rs, _ = event_horizon(M())
    @test (_rs .> 1.0) |> all
end

# test naked singularities
@test is_naked_singularity(KerrMetric(1.0, 0.0)) == false
@test is_naked_singularity(JohannsenPsaltisMetric(1.0, 0.998, 2.0)) == true
