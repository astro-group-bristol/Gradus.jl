using Test
using Gradus
using StaticArrays

# Tests to make sure the basic pointfunctions work

@testset "pointfunctions" begin
    u = @SVector [0.0, 100.0, deg2rad(85), 0.0]

    function run_pointfunction(m, pf)
        d = ThinDisc(10.0, 40.0)
        α, β, img = rendergeodesics(
            m,
            u,
            d,
            200.0,
            image_width = 100,
            image_height = 100,
            αlims = (-50, 50),
            βlims = (-50, 50),
            pf = pf,
            verbose = false,
        )
        img
    end

    @testset "redshift" begin
        # only implemented for the BoyerLindquist metrics at the moment
        for m in [KerrMetric(), KerrSpacetimeFirstOrder()]
            run_pointfunction(m, Gradus.ConstPointFunctions.redshift(m, u))
        end
        # smoke test passed
        @test true
    end
end
