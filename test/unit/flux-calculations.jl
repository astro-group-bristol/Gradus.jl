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

# proper area

# using Wikins+Fabian 2012 / Dauser+13 (though there is a typo in Dauser+13)
# even though their plot is correct
function area(a, r)
    A = r^4 + a^2 * r^2 + 2 * a^2 * r
    B = r^2 - 2 * r + a^2
    2 * π * √(A / B)
end

function proper_area(m, x)
    gcomp = Gradus.metric_components(m, x)
    det_g = √(gcomp[2] * gcomp[4])
    2π * det_g
end

area1 = map(rrange) do r
    x = SVector(r, π / 2)
    (2π * r) / proper_area(m, x)
end

area2 = map(rrange) do r
    (2π * r) / area(m.a, r)
end

@test area1 ≈ area2
