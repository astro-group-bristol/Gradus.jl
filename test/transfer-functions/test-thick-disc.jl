using Test
using Gradus

m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 10_000, deg2rad(75), 0.0)
d = ShakuraSunyaev(m)

tf = cunningham_transfer_function(m, x, d, 3.0; β₀ = 2.0)

total = sum(filter(!isnan, tf.f))
@test total ≈ 14.64279128586961 atol = 1e-4

m = KerrMetric(1.0, 0.2)
x = SVector(0.0, 10_000, deg2rad(20), 0.0)
d = ShakuraSunyaev(m; eddington_ratio = 0.2)

tf = cunningham_transfer_function(m, x, d, 5.469668466100368; β₀ = 2.0)
total = sum(filter(!isnan, tf.f))
@test total ≈ 21.581370829241525 atol = 1e-2

# the transfer function here is pretty horrible as it's almost impossible to actually
# see; this is a test to make sure it doesn't error
# an offset to the isco of 4-e2 resolves this, but that's quite a lot
tf = cunningham_transfer_function(m, x, d, Gradus.isco(m) + 1e-2; β₀ = 1.0)

function _th_ctx(;
    r_target,
    r = 10_000.0,
    th = 45,
    a = 0.998,
    edd_ratio = 0.3,
    β₀ = 1.0,
    kwargs...,
)
    m = KerrMetric(1.0, a)
    d = ShakuraSunyaev(m; eddington_ratio = edd_ratio)
    Gradus.cunningham_transfer_function(
        m,
        SVector(0.0, r, deg2rad(th), 0.0),
        d,
        r_target;
        β₀ = β₀,
        kwargs...,
    )
end

# various test cases
begin
    @time _th_ctx(; a = 0.0, r_target = 903.9954031222643, th = 70, β₀ = 1.5)
    @time _th_ctx(; a = 0.0, r_target = 6.0, th = 70, β₀ = 1.5)
    @time _th_ctx(; a = 0.0, r_target = 6.0, th = 45, β₀ = 1.0)
    @time _th_ctx(; a = 0.0, r_target = 6.0, th = 5, β₀ = 0.0)
end

# these have been problematic
begin
    _th_ctx(; r_target = 903.9954031222643, th = 85, β₀ = 1.5)
end
