using Test
using Gradus
using StaticArrays

# smoke test to make sure circular orbits work

@testset "circular-orbits" begin

    @testset "solve_equitorial_circular_orbit" begin
        # only implemented for the BoyerLindquist metrics at the moment
        # expected is sum of the circular orbit vϕ
        # last updated: 22 Apr 2022
        r_range = 6.0:0.5:10.0
        for (m, expected) in [
            (KerrMetric(M = 1.0, a = 0.0), 0.5432533297869712),
            (KerrMetric(M = 1.0, a = 1.0), 0.5016710246454921),
            (KerrMetric(M = 1.0, a = -1.0), 0.5993458160081419),
            (JohannsenMetric(M = 1.0, a = 1.0, α22 = 1.0), 0.4980454719932759),
        ]
            vϕs = solve_equitorial_circular_orbit(m, r_range)
            @test isapprox(sum(vϕs), expected, atol = 1e-6, rtol = 0.0)
        end
    end

    @testset "trace_equitorial_circular_orbit" begin
        m = KerrMetric(M = 1.0, a = 0.0)
        Gradus.isco(m)
        r_range = 6.0:0.1:10.0
        simsols = trace_equitorial_circular_orbit(m, r_range)
        # smoke test passed
        @test true
    end
end
