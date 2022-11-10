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

        intersected_simsols =
            filter(i -> i.prob.p.status == StatusCodes.IntersectedWithGeometry, simsols.u)
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
        # last computed 03/10/22: updated sampling method
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

    @testset "surrogate-disc-profiles" begin
        d = GeometricThinDisc(1.0, 40.0, deg2rad(90.0))

        x_train = SVector{2,Float64}[
            [-2.990522674851722, -0.758227051833861],
            [1.6957317742742792, 2.579567047821497],
            [0.49036042208502595, -3.04969081635049],
            [-2.42176062237317, 1.9203007984462452],
            [3.084878237683712, 0.2179644105725826],
            [-2.1448601698830245, -2.2034978296246726],
            [0.09483550185312209, 3.075782787915406],
            [2.007950237970996, -2.3347612107826397],
            [-3.0596834032923637, 0.3671029387860535],
            [2.5068396096880377, 1.7960603639659394],
        ]
        y_train = [
            0.21425436604784126,
            0.22201025331074273,
            0.22979387544026428,
            0.2377146955978877,
            0.2458179016185599,
            0.2532528211247887,
            0.2613657350968698,
            0.269506826643881,
            0.277668595159603,
            0.28600957447000247,
        ]
        # instantiate
        rbf = SurrogateDiscProfile(d, x_train, y_train)

        # test single evaluate
        @test all(evaluate(rbf, x_train[1]) .≈ y_train[1])
        # test many
        @test all(evaluate(rbf, x_train) .≈ y_train)
    end
end
