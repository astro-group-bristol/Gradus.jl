using Test
using Gradus
using StaticArrays

# Tests to make sure the basic different special radii functions work.

@testset "special-radii" begin
    all_metrics = (
        KerrSpacetimeFirstOrder(1.0, 0.998, 1.0),
        KerrSpacetimeFirstOrder(1.0, -0.998, 1.0),
        KerrMetric(1.0, 0.998),
        KerrMetric(1.0, -0.998),
        JohannsenMetric(M = 1.0, a = 0.998, α13 = 1.0),
        JohannsenMetric(M = 1.0, a = 0.998, α22 = 1.0),
        DilatonAxion(M = 1.0, a = 0.998, β = 0.2, b = 1.0),
    )

    @testset "iscos" begin
        # last computed 04/06/2022: tests with --math-mode=ieee
        for (m, expected) in zip(
            all_metrics,
            [
                1.2369706551751847,
                8.99437445480357,
                1.2369706551751847,
                8.99437445480357,
                2.8482863127671534,
                1.1306596884484472,
                6.880202293032178,
            ],
        )
            isco_r = Gradus.isco(m)
            @test isapprox(expected, isco_r; atol = 1e-5, rtol = 0.0)
        end
    end
end
