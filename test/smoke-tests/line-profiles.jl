
@testset "line-profiles" begin
    d = GeometricThinDisc(0.0, 300.0, π / 2)
    u = @SVector [0.0, 1000.0, deg2rad(40), 0.0]
    m = BoyerLindquistAD(M = 1.0, a = 0.998)

    @testset "cunningham" begin
        x = range(0.1, 1.2, 20)
        bins, lp = lineprofile(
            (re) -> re^(-3),
            m,
            u,
            d;
            num_points = 40,
            bins = x,
            num_re = 40,
            max_re = 50,
            finite_diff_order = 5,
        )
        @test isapprox(0.20282590591902194, sum(lp), atol = 1e-4)
    end
end
