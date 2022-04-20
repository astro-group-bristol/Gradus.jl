# Tests to make sure the basic `tracegeodesics` function works for (ideally) all metrics.
using Test, Gradus, StaticArrays

@testset "rendergeodesics" begin

    @testset "plain" begin
        u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
        for (m, expectation) in zip(
            [BoyerLindquistAD(), JohannsenAD(), BoyerLindquistFO(), MorrisThorneAD()],
            # expectation values for the sum of the image
            # last computed 20/04/2022: v0.1.0
            [8949.345629557964, 8949.36445058427, 8961.703066256105, 306.88361044534764],
        )
            img = rendergeodesics(
                m,
                u,
                200.0,
                fov_factor = 1.0,
                image_width = 100,
                image_height = 50,
                verbose = false,
            )
            image_fingerprint = sum(filter(!isnan, img))
            @test isapprox(expectation, image_fingerprint; atol=1.0, rtol=0.0)
        end
    end

    @testset "thin-disc" begin
        u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
        d = GeometricThinDisc(10.0, 40.0, deg2rad(90.0))
        for (m, expectation) in zip(
            [BoyerLindquistAD(), JohannsenAD(), BoyerLindquistFO(), MorrisThorneAD()],
            # expectation values for the sum of the image
            # last computed 20/04/2022: v0.1.0
            [86795.10968236077, 86796.73210488849, 81959.05505753662, 36803.48713530963],
        )
            img = rendergeodesics(
                m,
                u,
                d,
                200.0,
                fov_factor = 1.0,
                image_width = 100,
                image_height = 50,
                verbose = false,
            )
            image_fingerprint = sum(filter(!isnan, img)) 
            @test isapprox(expectation, image_fingerprint; atol=1.0, rtol=0.0)
        end
    end

end
