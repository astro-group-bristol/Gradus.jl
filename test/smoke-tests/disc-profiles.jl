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

        intersected_simsols =
            filter(i -> i.prob.p.status == StatusCodes.IntersectedWithGeometry, simsols.u)
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
        # last computed 03/10/22: updated sampling method
        expected_areas = [
            75.64302307389167,
            114.78829414215679,
            166.581361553731,
            263.3161757085493,
            700.3526293678786,
            739.2923798563061,
            645.4220326474972,
            477.8144794927597,
            359.1850275798407,
        ]
        @test all(isapprox.(areas1, expected_areas, atol = 1e-3))
        @test all(isapprox.(areas2, expected_areas, atol = 1e-3))
    end
end
