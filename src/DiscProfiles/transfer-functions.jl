function bin_transfer_function(time_delays, energy, flux; N = 300)
    energy_bins = range(extrema(energy)..., N)
    time_bins = range(extrema(time_delays)..., N)

    de = step(energy_bins)
    dt = step(time_bins)

    transfer_function = ThreadsX.map(Iterators.product(time_bins, energy_bins)) do (t, e)
        t_mask = @. (t ≤ time_delays) & (time_delays < (t + dt))
        e_mask = @. (e ≤ energy) & (energy < (e + de))
        sub_flux = @views(flux[t_mask.&e_mask])
        if length(sub_flux) > 0
            sum(sub_flux)
        else
            NaN
        end
    end

    time_bins, energy_bins, transfer_function
end
