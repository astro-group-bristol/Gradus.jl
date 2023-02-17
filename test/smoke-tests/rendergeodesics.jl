using Test
using Gradus
using StaticArrays

# Tests to make sure the basic `rendergeodesics` function works for (ideally) all metrics.

function _thick_disc(u)
    r = u[2]
    if r < 9.0 || r > 11.0
        return -1.0
    else
        x = r - 10.0
        sqrt(1 - x^2)
    end
end

@testset "plain" begin
    u = @SVector [0.0, 100.0, deg2rad(85), 0.0]

    # last computed 21/01/2023: shrink resolution
    expected = (8969.1564582409967, 8969.15634220181, 8977.502920124776, 413.4963434133726)
    result = map((
        KerrMetric(),
        JohannsenMetric(),
        KerrSpacetimeFirstOrder(),
        MorrisThorneWormhole(),
    )) do m
        α, β, img = rendergeodesics(
            m,
            u,
            200.0,
            fov_factor = 1.0,
            image_width = 20,
            image_height = 20,
            verbose = false,
        )
        image_fingerprint = sum(filter(!isnan, img))
    end
    for (e, v) in zip(expected, result)
        # this tolerance is kind of unacceptably high? todo: investigate why
        @test isapprox(e, v; rtol = 0.1)
    end
end

@testset "thin-disc" begin
    u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
    d = GeometricThinDisc(10.0, 40.0, deg2rad(90.0))

    # last computed 21/01/2023: shrink resolution
    expected = (29605.55590761622, 29605.556409711604, 29741.80749605271, 9858.77920909911)
    result = map((
        KerrMetric(),
        JohannsenMetric(),
        KerrSpacetimeFirstOrder(),
        MorrisThorneWormhole(),
    )) do m
        α, β, img = rendergeodesics(
            m,
            u,
            d,
            200.0,
            fov_factor = 1.0,
            image_width = 20,
            image_height = 20,
            verbose = false,
        )
        image_fingerprint = sum(filter(!isnan, img))
    end
    for (e, v) in zip(expected, result)
        # this tolerance is kind of unacceptably high? todo: investigate why
        @test isapprox(e, v; rtol = 0.1)
    end
end

@testset "shakura-sunyaev-disc" begin
    u = @SVector [0.0, 100.0, deg2rad(85), 0.0]

    # last computed 21/01/2023: shrink resolution
    expected = (34711.33445248479, 34711.33445255157)
    result = map((KerrMetric(), JohannsenMetric())) do m
        d = ShakuraSunyaev(m)
        α, β, img = rendergeodesics(
            m,
            u,
            d,
            200.0,
            fov_factor = 1.0,
            image_width = 20,
            image_height = 20,
            verbose = false,
        )
        image_fingerprint = sum(filter(!isnan, img))
    end
    for (e, v) in zip(expected, result)
        # this tolerance is kind of unacceptably high? todo: investigate why
        @test isapprox(e, v; rtol = 0.1)
    end
end

@testset "thick-disc" begin
    u = @SVector [0.0, 100.0, deg2rad(85), 0.0]
    d = ThickDisc(_thick_disc)

    # last computed 21/01/2023: shrink resolution
    expected = (16593.560393732, 16593.56001187974, 16847.84450997791, 5015.839213855068)
    result = map((
        KerrMetric(),
        JohannsenMetric(),
        KerrSpacetimeFirstOrder(),
        MorrisThorneWormhole(),
    )) do m
        α, β, img = rendergeodesics(
            m,
            u,
            d,
            200.0,
            fov_factor = 1.0,
            image_width = 20,
            image_height = 20,
            verbose = false,
        )
        image_fingerprint = sum(filter(!isnan, img))
    end
    for (e, v) in zip(expected, result)
        # this tolerance is kind of unacceptably high? todo: investigate why
        @test isapprox(e, v; rtol = 0.1)
    end
end
