using Test
using Gradus

function test_ctf1(a, angle, rₑ)
    m = KerrMetric(1.0, a)
    d = ThinDisc(0.0, Inf)
    x = SVector(0.0, 10_000, deg2rad(angle), 0.0)
    Gradus.cunningham_transfer_function(
        m,
        x,
        d,
        rₑ,
        ;
        chart = Gradus.chart_for_metric(m, 2 * x[2]),
        N = 80,
    )
end

function measure_ctf(ctf)
    sum(ctf.f .* ctf.g✶) / length(ctf.f)
end

# test for different angles
@test measure_ctf(test_ctf1(0.998, 3, 4.0)) ≈ 0.12078230790537134 atol = 1e-5
@test measure_ctf(test_ctf1(0.998, 35, 4.0)) ≈ 0.10309564736283175 atol = 1e-5
@test measure_ctf(test_ctf1(0.998, 74, 4.0)) ≈ 0.05346607396286273 atol = 1e-5
@test measure_ctf(test_ctf1(0.998, 85, 4.0)) ≈ 0.034591443363911956 atol = 1e-5

# different radii
@test measure_ctf(test_ctf1(0.998, 30, 4.0)) ≈ 0.10702748528655079 atol = 1e-5
@test measure_ctf(test_ctf1(0.998, 30, 7.0)) ≈ 0.11941254099997582 atol = 1e-5
@test measure_ctf(test_ctf1(0.998, 30, 10.0)) ≈ 0.12411523229602836 atol = 1e-5
@test measure_ctf(test_ctf1(0.998, 30, 15.0)) ≈ 0.1269951484304088 atol = 1e-5

# large radii
@test measure_ctf(test_ctf1(0.998, 30, 300.0)) ≈ 0.134116533963842 atol = 1e-5
@test measure_ctf(test_ctf1(0.998, 30, 800.0)) ≈ 0.1339736833006852 atol = 1e-5
@test measure_ctf(test_ctf1(0.998, 30, 1000.0)) ≈ 0.1359210533432173 atol = 1e-5

# ones that have been problematic in the past
# and should not raise any errors or warnings
test_ctf1(-0.6, 88, 784.8253509875607)
test_ctf1(-0.998, 88, 953.9915665264327)
test_ctf1(-0.450, 88, 952.1406350219423)
test_ctf1(0.0, 88, 952.1406350219423)
test_ctf1(0.0, 88, 631.1007589946363)
test_ctf1(0.9, 88, 952.1406350219423)
test_ctf1(0.0, 88, 950.8754196211696)
test_ctf1(0.744, 88, 3.1880132176627862)

# poor offset radius
test_ctf1(1.0, 88, 1.01)
