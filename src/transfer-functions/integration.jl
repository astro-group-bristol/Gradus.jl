g_to_g✶(g, gmin, gmax) = @. (g - gmin) / (gmax - gmin)
g✶_to_g(g✶, gmin, gmax) = @. (gmax - gmin) * g✶ + gmin

function interpolate_over_radii(itfs, g✶_grid)
    radii = map(itf -> itf.rₑ, itfs)

    if !issorted(radii)
        I = sortperm(radii)
        radii = radii[I]
        itfs = itfs[I]
    end

    fr_interp = map(g✶_grid) do g✶
        f = map(itfs) do itf
            itf.lower_f(g✶) + itf.upper_f(g✶)
        end
        DataInterpolations.LinearInterpolation(f, radii)
    end

    gmin_interp = DataInterpolations.LinearInterpolation(map(itf -> itf.gmin, itfs), radii)
    gmax_interp = DataInterpolations.LinearInterpolation(map(itf -> itf.gmax, itfs), radii)
    fr_interp, gmin_interp, gmax_interp
end

function integrate_edge(S, lim, lim_g✶, gmin, gmax)
    gh = g✶_to_g(lim_g✶, gmin, gmax)
    a = abs(√gh - √lim)

    2 * S(gh) * a
end

function _cunningham_integrand(f, g, gs, gmin, gmax)
    f * π * (g)^3 / (√(gs * (1 - gs)) * (gmax - gmin))
end

function integrate_bin(S, gmin, gmax, lo, hi; h = 2e-8)
    glo = clamp(lo, gmin, gmax)
    ghi = clamp(hi, gmin, gmax)

    lum = 0.0
    if glo == ghi
        return lum
    end

    g✶lo = g_to_g✶(lo, gmin, gmax)
    g✶hi = g_to_g✶(hi, gmin, gmax)

    if g✶lo < h
        if g✶hi > h
            lum += integrate_edge(S, glo, h, gmin, gmax)
            glo = g✶_to_g(h, gmin, gmax)
        else
            return integrate_edge(S, glo, g✶hi, gmin, gmax)
        end
    end

    if g✶hi > 1 - h
        if g✶lo < 1 - h
            lum += integrate_edge(S, ghi, 1 - h, gmin, gmax)
            ghi = g✶_to_g(1 - h, gmin, gmax)
        else
            return integrate_edge(S, ghi, g✶lo, gmin, gmax)
        end
    end

    res, _ = quadgk(S, glo, ghi)
    lum += res
    return lum
end

function _wrap_cunningham_interpolations(fr_interp, gmin, gmax, r, g✶_grid)
    f = map(fr_interp) do interp
        interp(r)
    end
    S = DataInterpolations.LinearInterpolation(f, g✶_grid)
    g -> begin
        g✶ = g_to_g✶(g, gmin, gmax)
        _cunningham_integrand(S(g✶), g, g✶, gmin, gmax)
    end
end

function _build_g✶_grid(Ng✶, h)
    g✶low = h
    g✶high = 1 - h
    map(1:Ng✶) do i
        g✶low + (g✶high - g✶low) / (Ng✶ - 1) * (i - 1)
    end
end

function integrate_drdg✶(
    ε,
    itfs::Vector{<:InterpolatedCunninghamTransferFunction{T}},
    radii,
    g_grid;
    Ng✶ = 10,
    h = 1e-8,
) where {T}
    # pre-allocate output
    flux = zeros(T, length(g_grid))

    # init params
    minrₑ, maxrₑ = extrema(radii)
    g✶_grid = _build_g✶_grid(Ng✶, h)
    fr_interp, gmin_interp, gmax_interp = interpolate_over_radii(itfs, g✶_grid)
    g_grid_view = @views g_grid[1:end-1]

    # build fine radial grid for trapezoidal integration
    fine_rₑ_grid = Grids._inverse_grid(minrₑ, maxrₑ, 1000) |> collect
    @inbounds for (i, rₑ) in enumerate(fine_rₑ_grid)
        gmin = gmin_interp(rₑ)
        gmax = gmax_interp(rₑ)
        S = _wrap_cunningham_interpolations(fr_interp, gmin, gmax, rₑ, g✶_grid)

        # trapezoidal integration weight
        if i == 1
            Δrₑ = (fine_rₑ_grid[i+1]) - (rₑ)
        elseif i == lastindex(fine_rₑ_grid)
            Δrₑ = (rₑ) - (fine_rₑ_grid[i-1])
        else
            Δrₑ = (fine_rₑ_grid[i+1]) - (fine_rₑ_grid[i-1])
        end

        weight = Δrₑ * rₑ * ε(rₑ)

        for j in eachindex(g_grid_view)
            glo = g_grid[j]
            ghi = g_grid[j+1]
            flux[j] += integrate_bin(S, gmin, gmax, glo, ghi; h = h) * weight
        end
    end

    Σflux = zero(T)
    @inbounds for i in eachindex(g_grid_view)
        glo = g_grid[i]
        ghi = g_grid[i+1]
        ḡ = (ghi + glo)
        flux[i] = flux[i] / ḡ
        Σflux += flux[i]
    end
    @. flux = flux / Σflux
    flux
end

export integrate_drdg✶, g✶_to_g, g_to_g✶
