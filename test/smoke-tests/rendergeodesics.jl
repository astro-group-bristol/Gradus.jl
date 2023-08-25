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

function _run_rendergeodesics(m, args...)
    x = SVector(0.0, 100.0, deg2rad(85), 0.0)
    α, β, img = rendergeodesics(
        m,
        x,
        args...,
        200.0;
        αlims = (-9.5, 9.5),
        βlims = (-9.5, 9.5),
        image_width = 20,
        image_height = 20,
        verbose = false,
    )
    image_fingerprint = sum(filter(!isnan, img))
end

metrics = (
    KerrMetric(),
    JohannsenMetric(),
    KerrSpacetimeFirstOrder(),
    MorrisThorneWormhole(),
    BumblebeeMetric(),
    KerrNewmanMetric(),
)

# shadows
# last computed 25/08/2023: change rendergeodesic axes behaviour
expected = [
    9009.452876609641,
    9009.448935932085,
    8873.083256336658,
    402.17907632733284,
    9009.452384885506,
    9009.451384824908,
]
result = [_run_rendergeodesics(m) for m in metrics]
# this tolerance is kind of unacceptably high? todo: investigate why
@test expected ≈ result rtol = 1e-1

# thin disc
d = GeometricThinDisc(0.0, 40.0, π / 2)
# last computed 25/08/2023: change rendergeodesic axes behaviour
expected = [
    38412.08347901267,
    38412.08386562321,
    37990.55152565702,
    9375.430228131403,
    38412.0832157869,
    38412.08517225652,
]
result = [_run_rendergeodesics(m, d) for m in metrics]
@test expected ≈ result rtol = 1e-1

# shakura sunyeav disc
result = map(metrics) do m
    if (m isa AbstractFirstOrderMetric) || (m isa MorrisThorneWormhole)
        0.0
    else
        _run_rendergeodesics(m, ShakuraSunyaev(m))
    end
end

result = [result...]
# last computed 25/08/2023: change rendergeodesic axes behaviour
expected =
    [34455.34416982827, 34455.344169980635, 0.0, 0.0, 34455.3441698318, 34455.34416971527]
@test expected ≈ result rtol = 1e-1

# thick disc
d = ThickDisc(_thick_disc)
result = [_run_rendergeodesics(m, d) for m in metrics]
# last computed 25/08/2023: change rendergeodesic axes behaviour
expected = [
    16918.69258396256,
    16918.689593279843,
    16920.323858977506,
    5104.032822512765,
    16918.692092917947,
    16918.691837255217,
]
@test expected ≈ result rtol = 1e-1
