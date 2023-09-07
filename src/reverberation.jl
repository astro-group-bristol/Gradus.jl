# first some utility functions

"""
Assumes x is equally spaced. Will not raise an error if not, but results
will be invalid.
"""
function extend_domain_with_zeros(x, y, x_max)
    Δx = x[2] - x[1]

    x̄ = range(minimum(x), x_max, step = Δx) |> collect
    ȳ = zeros(eltype(y), size(x̄))
    ȳ[firstindex(y):lastindex(y)] .= y

    x̄, ȳ
end

sum_impulse_response(f::AbstractMatrix) =
    sum(i -> isnan(i) ? zero(eltype(f)) : i, f, dims = 1) |> vec

phase_difference(𝔉ψ::AbstractArray{<:Complex}) = @. atan(imag(𝔉ψ) / (1 + real(𝔉ψ)))

function time_lag!(out, freq, ϕ)
    for i in eachindex(out)
        out[i] = ϕ[i] / (2π * freq[i])
    end
    out
end

function lag_frequency(t, ψ; R = 1)
    freq = FFTW.fftfreq(length(t), 1 / (t[2] - t[1]))
    𝔉ψ = R .* FFTW.fft(ψ)
    # pick only positive frequencies
    I = 1:(length(freq)÷2)
    ϕ = @views phase_difference(𝔉ψ[I])
    # store the time lag in ϕ directly
    time_lag!(ϕ, freq, ϕ)
    freq[I], -ϕ
end

function lag_frequency(t, f::AbstractMatrix; flo = 5e-5, kwargs...)
    ψ = sum_impulse_response(f)
    t_extended, ψ_extended = extend_domain_with_zeros(t, ψ, 1 / flo)

    lag_frequency(t_extended, ψ_extended; kwargs...)
end

function lag_frequency(
    m::AbstractMetric{T},
    x,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    Nr = 6000,
    bins = collect(range(0.0, 1.5, 500)),
    tbins = collect(range(0, 1000.0, 2000)),
    spectrum = PowerLawSpectrum(2),
    radii = collect(range(Gradus.isco(m) + 1e-2, 300.0, 100)),
    kwargs...,
) where {T}
    other_kwargs, em_setup = _EmissivityProfileSetup(T, spectrum; kwargs...)
    solver_kwargs, tf_setup = _TransferFunctionSetup(T; other_kwargs...)

    prof = emissivity_profile(em_setup, m, d, model; solver_kwargs...)
    t0 = continuum_time(m, x, model; solver_kwargs...)
    itb = interpolated_transfer_branches(tf_setup, m, x, d, radii; solver_kwargs...)

    flux = @time Gradus.integrate_lagtransfer(
        prof,
        itb,
        radii,
        bins,
        tbins;
        t0 = t0,
        Nr = Nr,
        h = tf_setup.h,
    )
    flux[flux.==0] .= NaN

    tbins, bins, flux
end

function continuum_time(m::AbstractMetric, x, model::AbstractCoronaModel; kwargs...)
    pos, _ = Gradus.sample_position_velocity(m, model)
    target = SVector(pos[2:end]...)
    _, _, gp, _ = Gradus.optimize_for_target(
        target,
        m,
        x;
        chart = Gradus.chart_for_metric(m, 2x[2]),
        callback = domain_upper_hemisphere(),
        kwargs...,
    )
    gp.x[1]
end

export lag_frequency, continuum_time
