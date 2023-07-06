using Test
using Gradus

# from Dauser+13, the Lorentz factor for the Kerr metric
# for Keplerian disc velocities in the equitorial plane
function _lorentz_factor(a, r)
    A = √(r^2 - 2r + a^2) * (r^(3 / 2) + a)
    B = √(r * √r + 2a - 3 * √r) * √(r^3 + a^2 * r + 2 * a^2) * r^(1 / 4)
    A / B
end

m = KerrMetric(1.0, 0.998)
d = GeometricThinDisc(0.0, 1000.0, π / 2)

rrange = collect(Gradus.Grids._geometric_grid(Gradus.isco(m), 1000.0, 100))

γ_gradus = map(rrange) do r
    x = SVector(0.0, r, π / 2, 0.0)
    Gradus.lorentz_factor(m, d, x)
end

γ_check = map(rrange) do r
    _lorentz_factor(m.a, r)
end

@test γ_gradus ≈ γ_check
