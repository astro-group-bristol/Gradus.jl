# Tests to make sure the basic pointfunctions work
using Test, Gradus, StaticArrays

@testset "pointfunctions" begin

    function run_pointfunction(m, pf)
        u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
        d = GeometricThinDisc(10.0, 40.0, deg2rad(90.0))
        img = rendergeodesics(
            m,
            u,
            d,
            200.0,
            fov_factor = 1.0,
            image_width = 100,
            image_height = 50,
            pf = pf,
            verbose = false,
        )
    end

    @testset "redshift" begin
        # only implemented for the BoyerLindquist metrics at the moment
        for m in [BoyerLindquistAD(), BoyerLindquistFO()]
            run_pointfunction(m, redshift)
        end
        # smoke test passed
        @test true
    end
end
