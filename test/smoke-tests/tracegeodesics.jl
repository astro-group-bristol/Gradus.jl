using Test
using Gradus
using StaticArrays

# Tests to make sure the basic `tracegeodesics` function works for (ideally) all metrics.

@testset "tracegeodesics" begin
    # tests if a single geodesic can be integrated
    function test_single(m, u, v)
        sol = tracegeodesics(m, u, v, (0.0, 200.0))
        @test sol.retcode == Gradus.SciMLBase.ReturnCode.Terminated
        sol
    end

    # tests if a multiple geodesics can be integrated
    function test_many(m, u, v)
        us = [u, u, u]
        vs = [v, v, v]
        simsols = tracegeodesics(m, us, vs, (0.0, 200.0))
        @test all(i -> simsols[i].retcode == Gradus.SciMLBase.ReturnCode.Terminated, 1:3)
        simsols
    end

    @testset "static-arrays" begin
        u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
        # arbitrary single velocity vector
        v = @SVector [0.0, -1.0, -0.0, -3.5e-6]
        for m in [
            KerrMetric(),
            JohannsenMetric(),
            KerrSpacetimeFirstOrder(),
            MorrisThorneWormhole(),
            BumblebeeMetric(),
            KerrNewmanMetric(),
        ]
            test_single(m, u, v)
            test_many(m, u, v)
            # smoke test passed
            @test true
        end
    end

    @testset "corona-models" begin
        # only implemented for KerrMetric at the moment
        # because of the vector to local sky methods
        m = KerrMetric(M = 1.0, a = 0.0)
        d = ThinDisc(Gradus.isco(m), 50.0)
        model = LampPostModel(h = 10.0, θ = deg2rad(0.001))

        for Sampler in [EvenSampler, WeierstrassSampler],
            Generator in [GoldenSpiralGenerator, RandomGenerator],
            Hemisphere in [LowerHemisphere, BothHemispheres]

            simsols = tracegeodesics(
                m,
                model,
                d,
                (0.0, 200.0);
                n_samples = 32,
                save_on = false,
                sampler = Sampler(generator = Generator(), domain = Hemisphere()),
            )
            # smoke test passed
            @test true
        end
    end
end
