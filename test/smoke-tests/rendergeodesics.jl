# Tests to make sure the basic `rendergeodesics` function works for (ideally) all metrics.

@testset "rendergeodesics" begin

    @testset "plain" begin
        u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
        for (m, expectation) in zip(
            [BoyerLindquistAD(), JohannsenAD(), BoyerLindquistFO(), MorrisThorneAD()],
            # expectation values for the sum of the image
            # last computed 04/06/2022: tests with --math-mode=ieee
            [8969.1564582409967, 8969.15634220181, 8977.502920124776, 413.49634341337264],
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
            # have to be really coarse cus the first order method is so variable???
            # the rest are very resolute
            @test isapprox(expectation, image_fingerprint; atol = 1e-3, rtol = 0.0)
        end
    end

    @testset "thin-disc" begin
        u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
        d = GeometricThinDisc(10.0, 40.0, deg2rad(90.0))
        for (m, expectation) in zip(
            [BoyerLindquistAD(), JohannsenAD(), BoyerLindquistFO(), MorrisThorneAD()],
            # expectation values for the sum of the image
            # last computed 04/06/2022: tests with --math-mode=ieee
            [86114.71322445248, 86233.47268379868, 81502.90270571726, 36043.47569730063],
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
            @test isapprox(expectation, image_fingerprint; atol = 1e-3, rtol = 0.0)
        end
    end

end
