
@testset "cunningham-transfer-functions" begin
    m = BoyerLindquistAD(1.0, 0.998)
    d = GeometricThinDisc(0.0, 400.0, deg2rad(90))

    # test for different angles
    for (angle, expected) in zip(
        [3, 35, 74, 85],
        # values last updated: 16/09/2022 - reduced number of points in test set
        [49.09620827339265, 45.16844822761859, 27.107307569529434, 19.664521288378978],
    )
        u = @SVector [0.0, 1000.0, deg2rad(angle), 0.0]
        ctf = cunningham_transfer_function(
            m,
            u,
            d,
            4.0,
            2000.0;
            finite_diff_order = 5,
            num_points = 200,
        )
        s = sum(ctf.f)
        @test isapprox(s, expected, atol = 1e-2, rtol = 0.0)
    end

    # different radii
    for (r, expected) in zip(
        [4.0, 7.0, 10.0, 15.0],
        # values last updated: 16/09/2022 - reduced number of points in test set
        [46.446844434569144, 51.08146707424876, 52.788771093295686, 53.78540573544522],
    )
        u = @SVector [0.0, 1000.0, deg2rad(30), 0.0]
        ctf = cunningham_transfer_function(
            m,
            u,
            d,
            r,
            2000.0;
            finite_diff_order = 5,
            num_points = 200,
        )
        s = sum(ctf.f)
        @test isapprox(s, expected, atol = 1e-2, rtol = 0.0)
    end
end
