# Tests to make sure the basic different special radii functions work.

@testset "special-radii" begin
    all_metrics = (
        BoyerLindquistFO(1.0, 0.998, 1.0),
        BoyerLindquistFO(1.0, -0.998, 1.0),
        BoyerLindquistAD(1.0, 0.998),
        BoyerLindquistAD(1.0, -0.998),
        JohannsenAD(M = 1.0, a = 0.998, α13 = 1.0),
        JohannsenAD(M = 1.0, a = 0.998, α22 = 1.0),
        DilatonAxionAD(M = 1.0, a = 0.998, β = 0.2, b = 1.0),
    )

    @testset "iscos" begin
        for m in all_metrics
            @test 100 > Gradus.isco(m) > 0.1
        end
    end
end
