using Test
using Gradus
using StaticArrays

# smoke test to make sure vornoi tesselation works

@testset "disc-profiles" begin
    @testset "voronoi-tesselation" begin
        m = KerrMetric(M = 1.0, a = 1.0)
        d = ThinDisc(1.0, 40.0)
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
        areas1 = Gradus.getareas(vdp1)
        areas2 = Gradus.getareas(vdp2)

        # values computed under visual inspection
        # last computed 21/06/23: updated sampling method
        expected_areas = [
            75.85149246031585,
            115.17768867370118,
            167.42061763426906,
            266.0208303714845,
            700.9761598342642,
            740.9728444986962,
            644.1965787516679,
            476.3536856817637,
            355.53014995977816,
        ]

        @test areas1 ≈ expected_areas atol = 1e-3
        @test areas2 ≈ expected_areas atol = 1e-3
    end
end
