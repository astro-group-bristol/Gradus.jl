using Test, Gradus

function _run_prerender(m, args...)
    x = SVector(0.0, 100.0, deg2rad(85), 0.0)
    α, β, cache = prerendergeodesics(
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
    pf =
        PointFunction((m, gp, λ_max) -> gp.λ_max) ∘
        FilterPointFunction((m, gp, λ_max) -> gp.λ_max < λ_max, NaN)
    img = apply(pf, cache)
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

# last computed 25/08/2023: change rendergeodesic axes behaviour
expected = [
    9009.452876609641,
    9009.448935932085,
    8873.083256336658,
    402.17907632733284,
    9009.452384885506,
    9009.451384824908,
]
result = [_run_prerender(m) for m in metrics]
@test expected ≈ result rtol = 1e-1
