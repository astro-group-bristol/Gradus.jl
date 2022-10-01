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

    transfer_function = ThreadsX.map(Iterators.product(time_bins, energy_bins)) do (t, e)
        t_mask = @. (t ≤ time_delays) & (time_delays < (t + dt))
        e_mask = @. (e ≤ energy) & (energy < (e + de))
        sub_flux = @views(flux[t_mask.&e_mask])
        if length(sub_flux) > 0
            sum(sub_flux) / (dt * de)
        else
            NaN
        end
    end

    time_bins, energy_bins, transfer_function
end

export bin_transfer_function
