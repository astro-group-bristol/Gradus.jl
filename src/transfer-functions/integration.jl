g_to_g✶(g, gmin, gmax) = @. (g - gmin) / (gmax - gmin)
g✶_to_g(g✶, gmin, gmax) = @. (gmax - gmin) * g✶ + gmin

function integrate_edge(S, lim, lim_g✶, gmin, gmax)
    gh = g✶_to_g(lim_g✶, gmin, gmax)
    a = abs(√gh - √lim)
    2 * S(gh) * a
end

@fastmath function _cunningham_integrand(f, g, gs, gmin, gmax)
    f * π * (g)^3 / (√(gs * (1 - gs)) * (gmax - gmin))
end

function integrate_bin(S, gmin, gmax, lo, hi; h = 2e-8, segbuf = nothing)
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

    res, _ = quadgk(S, glo, ghi; segbuf = segbuf)
    lum += res
    return lum
end

function _transform_g✶_to_g(branch)
    gmin = branch.gmin
    gmax = branch.gmax

    function _g_lower_f(g)
        g✶ = g_to_g✶(g, gmin, gmax)
        _cunningham_integrand(branch.lower_f(g✶), g, g✶, gmin, gmax)
    end
    function _g_upper_f(g)
        g✶ = g_to_g✶(g, gmin, gmax)
        _cunningham_integrand(branch.upper_f(g✶), g, g✶, gmin, gmax)
    end
    function _g_lower_t(g)
        branch.lower_t(g_to_g✶(g, gmin, gmax))
    end
    function _g_upper_t(g)
        branch.upper_t(g_to_g✶(g, gmin, gmax))
    end

    _g_lower_f, _g_upper_f, _g_lower_t, _g_upper_t
end

function _build_g✶_grid(Ng✶, h)
    g✶low = h
    g✶high = 1 - h
    map(1:Ng✶) do i
        g✶low + (g✶high - g✶low) / (Ng✶ - 1) * (i - 1)
    end
end

function _trapezoidal_weight(X, x, i)
    if i == 1
        X[i+1] - x
    elseif i == lastindex(X)
        x - X[i-1]
    else
        X[i+1] - X[i-1]
    end
end

function _normalize(flux::AbstractArray{T}, grid) where {T}
    Σflux = zero(T)
    @inbounds for i = 1:length(grid)-1
        glo = grid[i]
        ghi = grid[i+1]
        ḡ = (ghi + glo)
        flux[i] = flux[i] / ḡ
        Σflux += flux[i]
    end
    @. flux = flux / Σflux
    flux
end

function integrate_drdg(
    ε,
    itb::InterpolatingTransferBranches{T},
    radii,
    g_grid;
    h = 1e-8,
) where {T}
    # pre-allocate output
    flux = zeros(T, length(g_grid))

    # init params
    minrₑ, maxrₑ = extrema(radii)
    g_grid_view = @views g_grid[1:end-1]

    # build fine radial grid for trapezoidal integration
    fine_rₑ_grid = Grids._inverse_grid(minrₑ, maxrₑ, 1000) |> collect

    # allocate a segbuf for Gauss-Kronrod
    segbuf = alloc_segbuf(Float64)

    @inbounds for (i, rₑ) in enumerate(fine_rₑ_grid)
        branch = itb(rₑ)
        F1, F2, T1, T2 = _transform_g✶_to_g(branch)

        Δrₑ = _trapezoidal_weight(fine_rₑ_grid, rₑ, i)
        weight = Δrₑ * rₑ * ε(rₑ)

        for j in eachindex(g_grid_view)
            glo = g_grid[j]
            ghi = g_grid[j+1]
            f1 =
                integrate_bin(
                    F1,
                    branch.gmin,
                    branch.gmax,
                    glo,
                    ghi;
                    h = h,
                    segbuf = segbuf,
                ) * weight
            f2 =
                integrate_bin(
                    F2,
                    branch.gmin,
                    branch.gmax,
                    glo,
                    ghi;
                    h = h,
                    segbuf = segbuf,
                ) * weight

            flux[j] += f1 + f2
        end
    end

    _normalize(flux, g_grid)
end

export integrate_drdg, g✶_to_g, g_to_g✶
