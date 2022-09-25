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
