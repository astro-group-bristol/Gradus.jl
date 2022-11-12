
@testset "line-profiles" begin
    d = GeometricThinDisc(0.0, 300.0, π / 2)
    u = @SVector [0.0, 1000.0, deg2rad(40), 0.0]
    m = BoyerLindquistAD(M = 1.0, a = 0.998)

    @testset "cunningham" begin
        x = range(0.1, 1.2, 20)
        bins = range(0.0, 1.5, 200)
        _, lp = lineprofile(bins, (r) -> r^(-3), m, u, d; N = 40, numrₑ = 30)
        @test isapprox(1.0, sum(lp), atol = 1e-4)
    end
end
