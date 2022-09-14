# smoke test to make sure vornoi tesselation works

@testset "disc-profiles" begin

    @testset "voronoi-tesselation" begin
        m = BoyerLindquistAD(M = 1.0, a = 1.0)
        d = GeometricThinDisc(1.0, 40.0, deg2rad(90.0))
        model = LampPostModel(h = 10.0, θ = deg2rad(0.0001))
        simsols = tracegeodesics(
            m,
            model,
            d,
            (0.0, 200.0),
            n_samples = 10,
            sampler = EvenSampler(domain = LowerHemisphere()),
        )

        intersected_simsols = filter(i -> i.retcode == :Intersected, simsols.u)
        sd_endpoints = map(sol -> getgeodesicpoint(m, sol), intersected_simsols)

        # test ensemble solution constructor
        vdp1 = VoronoiDiscProfile(m, d, simsols)
        @test true

        # test endpoint constructor
        vdp2 = VoronoiDiscProfile(m, d, sd_endpoints)
        @test true

        # check areas
        areas1 = getareas(vdp1)
        areas2 = getareas(vdp2)

        # values computed under visual inspection
        # last computed 15/11/2022: continuous callback for disc intersection
        expected_areas = [
            89.79326183804586,
            117.69063502043767,
            141.71137385746914,
            233.49741646976452,
            683.6491180876137,
            736.059111143716,
            675.7846796957201,
            490.44994928380703,
            372.2875995001252,
        ]
        @test all(isapprox.(areas1, expected_areas, atol = 1e-3))
        @test all(isapprox.(areas2, expected_areas, atol = 1e-3))

    end
end
