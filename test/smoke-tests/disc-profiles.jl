using Test
using Gradus
using StaticArrays

# smoke test to make sure vornoi tesselation works

@testset "disc-profiles" begin
    @testset "voronoi-tesselation" begin
        m = KerrMetric(M = 1.0, a = 1.0)
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

        intersected_simsols = filter(
            i ->
                Gradus.get_status_code(i.prob.p) == StatusCodes.IntersectedWithGeometry,
            simsols.u,
        )
        sd_endpoints = map(sol -> unpack_solution(m, sol), intersected_simsols)

        # test ensemble solution constructor
        vdp1 = VoronoiDiscProfile(m, d, simsols)
        @test true

        # test endpoint constructor
        vdp2 = VoronoiDiscProfile(m, d, sd_endpoints)
        @test true

        # check areas
        areas1 = getareas(vdp1)
        areas2 = getareas(vdp2)

        # values computed under visual inspection
        # last computed 21/06/23: updated sampling method
        expected_areas = [
            75.67457406264847,
            114.75628423611069,
            166.48330049359774,
            263.58121169091913,
            700.3117449256706,
            739.4180350022357,
            645.1747233625393,
            478.1614271222751,
            358.7139818627764,
        ]

        @test areas1 ≈ expected_areas atol = 1e-3
        @test areas2 ≈ expected_areas atol = 1e-3
    end
end
