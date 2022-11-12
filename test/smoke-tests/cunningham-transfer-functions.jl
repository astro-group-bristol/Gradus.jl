
@testset "cunningham-transfer-functions" begin
    m = BoyerLindquistAD(1.0, 0.998)
    d = GeometricThinDisc(0.0, 400.0, deg2rad(90))

    # test for different angles
    for (angle, expected) in zip(
        [3, 35, 74, 85],
        # values last updated: 10/11/2022 - reimplemented algorithm
        [24.953515365137005, 22.91076063485734, 13.795031677382394, 11.086660268049126],
    )
        u = @SVector [0.0, 1000.0, deg2rad(angle), 0.0]
        ctf = cunningham_transfer_function(m, u, d, 4.0)
        s = sum(ctf.f)
        @test isapprox(s, expected, atol = 1e-1, rtol = 0.0)
    end

    # different radii
    for (r, expected) in zip(
        [4.0, 7.0, 10.0, 15.0],
        # values last updated: 10/11/2022 - reimplemented algorithm
        [23.491668065237214, 25.997534896876644, 26.73859117429336, 27.186165046726007],
    )
        u = @SVector [0.0, 1000.0, deg2rad(30), 0.0]
        ctf = cunningham_transfer_function(m, u, d, r)
        s = sum(ctf.f)
        @test isapprox(s, expected, atol = 1e-2, rtol = 0.0)
    end
end
