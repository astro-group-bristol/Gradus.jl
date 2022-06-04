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
        # last computed: 29/04/2022
        expected_areas = [
            92.71878363579425,
            122.84126831622028,
            150.79247196704486,
            251.36539040790862,
            690.909009896099,
            749.6987472036542,
            666.1346645350133,
            472.7915541962405,
            352.20220398142754,
        ]

        @show areas1
        @show areas2

        @test all(isapprox.(areas1, expected_areas, atol = 1e-6))
        @test all(isapprox.(areas2, expected_areas, atol = 1e-6))

    end
end
