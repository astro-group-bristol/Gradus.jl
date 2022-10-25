# Tests to make sure the basic `rendergeodesics` function works for (ideally) all metrics.

function _thick_disc(u)
    r = u[2]
    if r < 9.0 || r > 11.0
        return -1.0
    else
        x = r - 10.0
        sqrt(1 - x^2)
    end
end

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
            @test isapprox(expectation, image_fingerprint; atol = 5.0, rtol = 0.0)
        end
    end

    @testset "thin-disc" begin
        u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
        d = GeometricThinDisc(10.0, 40.0, deg2rad(90.0))
        for (m, expectation) in zip(
            [BoyerLindquistAD(), JohannsenAD(), BoyerLindquistFO(), MorrisThorneAD()],
            # expectation values for the sum of the image
            # last computed 15/11/2022: continuous callback for disc intersection
            [90629.67068537966, 90558.28684240118;, 86277.8210033628, 39094.80412600604],
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
            # this tolerance is kind of unacceptably high? todo: investigate why
            @test isapprox(expectation, image_fingerprint; atol = 5.0, rtol = 0.0)
        end
    end

    @testset "shakura-sunyaev-disc" begin
        u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
        for (m, expectation) in zip(
            [BoyerLindquistAD(), JohannsenAD()],
            # expectation values for the sum of the image
            # last computed 11/10/2022: initial values
            [282541.415414396, 282541.4154167309],
        )
            d = ShakuraSunyaev(m)
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
            @test isapprox(expectation, image_fingerprint; atol = 5.0, rtol = 0.0)
        end
    end

    @time @testset "thick-disc" begin
        u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
        d = ThickDisc(_thick_disc)
        for (m, expectation) in zip(
            [BoyerLindquistAD(), JohannsenAD(), BoyerLindquistFO(), MorrisThorneAD()],
            # expectation values for the sum of the image
            # last computed 11/10/2022: initial values
            [17755.784049024965, 17755.7836819751439, 18018.50638877236, 5401.369964242934],
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
            @test isapprox(expectation, image_fingerprint; atol = 5.0, rtol = 0.0)
        end
    end

end
