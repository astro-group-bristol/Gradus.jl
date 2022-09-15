
@testset "cunningham-transfer-functions" begin
    m = BoyerLindquistAD(1.0, 0.998)
    d = GeometricThinDisc(0.0, 400.0, deg2rad(90))

    # test for different angles
    for (angle, expected) in zip(
        [3, 35, 74, 85],
        [248.95547957677172, 228.40998327542954, 136.93067378322644, 99.95201562448777],
    )
        u = @SVector [0.0, 1000.0, deg2rad(angle), 0.0]
        ctf = Gradus.cunningham_transfer_function(
            m,
            u,
            d,
            4.0,
            2000.0;
            finite_diff_order = 5,
        )
        s = sum(ctf.f)
        @test isapprox(s, expected, atol = 1e-2, rtol = 0.0)
    end

    # different radii
    for (r, expected) in zip(
        [4.0, 7.0, 10.0, 15.0],
        [234.24090499819158, 259.3365885895835, 266.8322238866137, 271.2289331188517],
    )
        u = @SVector [0.0, 1000.0, deg2rad(30), 0.0]
        ctf = Gradus.cunningham_transfer_function(
            m,
            u,
            d,
            r,
            2000.0;
            finite_diff_order = 5,
        )
        s = sum(ctf.f)
        @test isapprox(s, expected, atol = 1e-2, rtol = 0.0)
    end
end
