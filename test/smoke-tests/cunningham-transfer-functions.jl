using Test
using Gradus
using StaticArrays

@testset "cunningham-transfer-functions" begin
    m = KerrMetric(1.0, 0.998)
    d = GeometricThinDisc(0.0, 400.0, deg2rad(90))

    # test for different angles
    for (angle, expected) in zip(
        [3, 35, 74, 85],
        # values last updated: 24/04/2023
        [
            0.24325323183855346,
            0.21844565584667297,
            0.11478208978397744,
            0.06742095227866991,
        ],
    )
        u = @SVector [0.0, 1000.0, deg2rad(angle), 0.0]
        ctf = cunningham_transfer_function(m, u, d, 4.0)
        s = sum(ctf.f) / length(ctf.f)
        @test isapprox(s, expected, atol = 1e-3, rtol = 0.0)
    end

    # different radii
    for (r, expected) in zip(
        [4.0, 7.0, 10.0, 15.0],
        # values last updated: 24/04/2023
        [0.22097479890817756, 0.2599801835844643, 0.26738568957410025, 0.2718492609857089],
    )
        u = @SVector [0.0, 1000.0, deg2rad(30), 0.0]
        ctf = cunningham_transfer_function(m, u, d, r)
        s = sum(ctf.f) / length(ctf.f)
        @test isapprox(s, expected, atol = 1e-2, rtol = 0.0)
    end
end
