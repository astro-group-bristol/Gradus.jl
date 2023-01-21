using Test
using Gradus
using StaticArrays

@testset "cunningham-transfer-functions" begin
    m = KerrSpacetime(1.0, 0.998)
    d = GeometricThinDisc(0.0, 400.0, deg2rad(90))

    # test for different angles
    for (angle, expected) in zip(
        [3, 35, 74, 85],
        # values last updated: 12/11/2022 - use means
        [0.2495227986998514, 0.22911943722172406, 0.13792349421541406, 0.10360637290617522],
    )
        u = @SVector [0.0, 1000.0, deg2rad(angle), 0.0]
        ctf = cunningham_transfer_function(m, u, d, 4.0)
        s = sum(ctf.f) / length(ctf.f)
        @test isapprox(s, expected, atol = 1e-3, rtol = 0.0)
    end

    # different radii
    for (r, expected) in zip(
        [4.0, 7.0, 10.0, 15.0],
        # values last updated: 12/11/2022 - use means
        [0.23492490558931373, 0.2599801835844643, 0.26738568957410025, 0.2718492609857089],
    )
        u = @SVector [0.0, 1000.0, deg2rad(30), 0.0]
        ctf = cunningham_transfer_function(m, u, d, r)
        s = sum(ctf.f) / length(ctf.f)
        @test isapprox(s, expected, atol = 1e-2, rtol = 0.0)
    end
end
