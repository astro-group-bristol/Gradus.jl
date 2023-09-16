using Test
using Gradus

function calculate_2d_transfer_function(m, x, model, itb, prof, radii)
    bins = collect(range(0.0, 1.5, 100))
    tbins = collect(range(0, 100.0, 100))

    t0 = continuum_time(m, x, model)

    flux = Gradus.integrate_lagtransfer(
        prof,
        itb,
        radii,
        bins,
        tbins;
        t0 = t0,
        Nr = 100,
        h = 1e-8,
    )

    flux[flux.==0] .= NaN
    bins, tbins, flux
end

function calc_lag_freq(m, d, model, radii, itb)
    prof = emissivity_profile(m, d, model; n_samples = 500)
    E, t, f = calculate_2d_transfer_function(m, x, model, itb, prof, radii)
    freq, τ = lag_frequency(t, f)
end

m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 10_000.0, deg2rad(45), 0.0)
d = ThinDisc(0.0, Inf)

radii = Gradus.Grids._inverse_grid(Gradus.isco(m), 100.0, 10)
itb = Gradus.interpolated_transfer_branches(m, x, d, radii; β₀ = 2.0)

model1 = LampPostModel()
freq1, τ1 = calc_lag_freq(m, d, model1, radii, itb)

# smoke test
@test sum(freq1) ≈ 2449.8787687490535 atol = 1e-2
@test τ1[132] ≈ 9.281677930459137 atol = 1e-2

# check it works for thick discs too
d_thick = ShakuraSunyaev(m)
itb_thick = Gradus.interpolated_transfer_branches(m, x, d_thick, radii; β₀ = 2.0)

freq2, τ2 = calc_lag_freq(m, d_thick, model1, radii, itb_thick)

# should be the same at this inclination
@test sum(freq2) ≈ 2449.8787687490535 atol = 1e-2
