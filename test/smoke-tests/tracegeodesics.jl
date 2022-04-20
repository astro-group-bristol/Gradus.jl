# Tests to make sure the basic `tracegeodesics` function works for (ideally) all metrics.
using Test, Gradus, StaticArrays

@testset "tracegeodesics" begin
    # tests if a single geodesic can be integrated
    function test_single(m, u, v)
        sol = tracegeodesics(m, u, v, (0.0, 200.0))
        @test sol.retcode == :Terminated
        sol
    end

    # tests if a multiple geodesics can be integrated
    function test_many(m, u, v)
        us = [u, u, u]
        vs = [v, v, v]
        simsols = tracegeodesics(m, us, vs, (0.0, 200.0))
        @test all(i -> simsols[i].retcode == :Terminated, 1:3)
        simsols
    end

    @testset "static-arrays" begin
        u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
        # arbitrary single velocity vector
        v = @SVector [0.0, -1.0, -0.0, -3.5e-6]
        for m in [BoyerLindquistAD(), JohannsenAD(), BoyerLindquistFO(), MorrisThorneAD()]
            test_single(m, u, v)
            test_many(m, u, v)
        end
    end

    @testset "regular-arrays" begin
        u = Float64[0.0, 100.0, deg2rad(85), 0.0]
        # arbitrary single velocity vector
        v = Float64[0.0, -1.0, -0.0, -3.5e-6]
        for m in [BoyerLindquistAD(), JohannsenAD(), BoyerLindquistFO(), MorrisThorneAD()]
            test_single(m, u, v)
            test_many(m, u, v)
        end
    end
end
