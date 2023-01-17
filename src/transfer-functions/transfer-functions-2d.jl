function bin_transfer_function(
    time_delays,
    energy,
    flux;
    N_E = 300,
    N_t = 300,
    energy_lims = extrema(energy),
    time_lims = extrema(time_delays),
)
    energy_bins = range(energy_lims..., N_E)
    time_bins = range(time_lims..., N_t)

    de = step(energy_bins)
    dt = step(time_bins)

    transfer_function =
        bucket(energy, time_delays, flux, energy_bins, time_bins; reduction = sum)

    @. transfer_function = transfer_function / (de * dt)
    transfer_function[transfer_function.==0.0] .= NaN

    time_bins, energy_bins, transfer_function
end

function radial_sampler(m::AbstractMetricParams, u, Nr::Int, Nφ::Int; rmax = 400, φmax = 2π)
    rs = Grids.geometric_grid(1.0, rmax, Nr)
    φs = range(0.0, φmax, Nφ)

    αs = vec([r * cos(φ) for r in rs, φ in φs])
    βs = vec([r * sin(φ) * sin(u[3]) for r in rs, φ in φs])

    f = (i) -> begin
        α = αs[i]
        β = βs[i]
        map_impact_parameters(m, u, α, β)
    end
    αs, βs, f
end

function bin_and_interpolate(
    X,
    y::AbstractArray{T};
    log_bins = false,
    nbins = 1000,
    reduction = mean,
) where {T}
    bins = if log_bins
        10 .^ range(log10(minimum(X)), log10(maximum(X)), nbins)
    else
        range(minimum(X), maximum(X), nbins)
    end

    y_binned = bucket(X, y, bins; reduction = reduction)
    DataInterpolations.LinearInterpolation(y_binned, bins)
end

export bin_transfer_function
