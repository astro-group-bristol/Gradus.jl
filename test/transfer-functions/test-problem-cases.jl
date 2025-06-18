using Test, Gradus

# Utility functions
_x(r, th) = SVector(0.0, r, deg2rad(th), 0.0)
_ctx(; a, r, th, r_target, d = ThinDisc(0.0, Inf), kwargs...) =
    Gradus.cunningham_transfer_function(
        KerrMetric(1.0, a),
        _x(r, th),
        d,
        r_target;
        kwargs...,
    )
_fulltf(; a, r, th, d = ThinDisc(0.0, Inf), kwargs...) =
    Gradus.transferfunctions(KerrMetric(1.0, a), _x(r, th), d; kwargs...)

# These are transfer functions that have caused issue in the past
_ctx(; a = 0.998, r = 500_000.0, th = 88.0, r_target = 1.2469706551751847)
_ctx(;
    a = 0.10324137931034483,
    r = 500_000.0,
    th = 82.06896551724138,
    r_target = 21.755193176415617,
)
_ctx(; a = 0.0, r = 500_000.0, th = 88.0, r_target = 264.549754423346)
_ctx(; a = 0.998, r = 500_000.0, th = 88.0, r_target = 1.2369706551751847)
_ctx(; a = 0.034413793103448276, r = 500_000.0, th = 88.0, r_target = 396.93135746662)
_ctx(; a = 0.034413793103448276, r = 500_000.0, th = 88.0, r_target = 377.0698611)
_ctx(; a = 0.034413793103448276, r = 500_000.0, th = 88.0, r_target = 417.83902340237086)
_ctx(; a = 0.0, r = 500_000.0, th = 88.0, r_target = 794.4185036834359)
_ctx(; a = 0.9291724137931034, r = 500_000.0, th = 88.0, r_target = 2.1204839212537308)

# _fulltf(; a = 0.10324137931034483, r = 500_000.0, th = 82.06896551724138, verbose = true)

# _fulltf(; a = 0.034413793103448276, r = 500_000.0, th = 88.0, verbose = true)
