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

function _transform_g✶_to_g(branch, h)
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
        g✶ = max(h, min(1 - h, g_to_g✶(g, gmin, gmax)))
        branch.lower_t(g✶)
    end
    function _g_upper_t(g)
        g✶ = max(h, min(1 - h, g_to_g✶(g, gmin, gmax)))
        branch.upper_t(g✶)
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

function _normalize(flux::AbstractVector{T}, grid) where {T}
    Σflux = zero(T)
    @inbounds for i = 1:length(grid)-1
        ḡ = (grid[i+1] + grid[i])
        flux[i] = flux[i] / ḡ
        Σflux += flux[i]
    end
    @. flux = flux / Σflux
    flux
end

function _normalize(flux::AbstractMatrix{T}, grid) where {T}
    Σflux = zero(T)
    @views @inbounds for i = 1:length(grid)-1
        ḡ = (grid[i+1] + grid[i])
        @. flux[i,:] = flux[i,:] / ḡ
        Σflux += sum(flux[i, :])
    end
    @. flux = flux / Σflux
    flux
end

function integrate_lineprofile(
    prof::AbstractDiscProfile,
    itb::InterpolatingTransferBranches,
    radii,
    g_grid;
    kwargs...,
)
    integrate_lineprofile(prof.f.ε, itb, radii, g_grid; kwargs...)
end

function integrate_lineprofile(
    ε,
    itb::InterpolatingTransferBranches{T},
    radii,
    g_grid;
    Nr = 1000,
    kwargs...,
) where {T}
    # pre-allocate output
    flux = zeros(T, length(g_grid))
    fine_rₑ_grid = Grids._inverse_grid(extrema(radii)..., Nr) |> collect
    _integrate_fluxbin!(flux, ε, itb, fine_rₑ_grid, g_grid; kwargs...)
    _normalize(flux, g_grid)
end

function integrate_lagtransfer(
    profile::AbstractDiscProfile,
    itb::InterpolatingTransferBranches{T},
    radii,
    g_grid,
    t_grid;
    Nr = 1000,
    kwargs...,
) where {T}
    # pre-allocate output
    flux = zeros(T, (length(g_grid), length(t_grid)))
    fine_rₑ_grid = Grids._geometric_grid(extrema(radii)..., Nr) |> collect
    _integrate_fluxbin!(flux, profile, itb, fine_rₑ_grid, g_grid, t_grid; kwargs...)
    _normalize(flux, g_grid)
end

function _integrate_fluxbin!(
    flux::AbstractVector,
    ε,
    itb::InterpolatingTransferBranches{T},
    radii_grid,
    g_grid;
    h = 1e-8,
    # allocate a segbuf for Gauss-Kronrod
    segbuf = alloc_segbuf(T),
) where {T}
    g_grid_view = @views g_grid[1:end-1]
    # build fine radial grid for trapezoidal integration
    @inbounds for (i, rₑ) in enumerate(radii_grid)
        branch = itb(rₑ)
        F1, F2, _, _ = _transform_g✶_to_g(branch, h)
        Δrₑ = _trapezoidal_weight(radii_grid, rₑ, i)
        # integration weight for this annuli
        θ = Δrₑ * rₑ * ε(rₑ)
        @inbounds for j in eachindex(g_grid_view)
            glo = g_grid[j]
            ghi = g_grid[j+1]
            f = integrate_bin(
                g -> F1(g) + F2(g),
                branch.gmin,
                branch.gmax,
                glo,
                ghi;
                h = h,
                segbuf = segbuf,
            )
            flux[j] += f * θ
        end
    end
end

function _integrate_fluxbin!(
    flux::AbstractMatrix,
    profile,
    itb::InterpolatingTransferBranches{T},
    radii_grid,
    g_grid,
    t_grid;
    h = 1e-8,
    # allocate a segbuf for Gauss-Kronrod
    segbuf = alloc_segbuf(T),
    t0 = 0,
) where {T}
    g_grid_view = @views g_grid[1:end-1]
    # build fine radial grid for trapezoidal integration
    @inbounds for (i, rₑ) in enumerate(radii_grid)
        branch = itb(rₑ)
        F1, F2, T1, T2 = _transform_g✶_to_g(branch, h)
        Δrₑ = _trapezoidal_weight(radii_grid, rₑ, i)

        # integration weight for this annuli
        θ = Δrₑ * rₑ * profile.f.ε(rₑ)
        # time delay for this annuli
        t_source_disc = profile.t.t(rₑ)

        @inbounds for j in eachindex(g_grid_view)
            glo = g_grid[j]
            ghi = g_grid[j+1]
            f1 = integrate_bin(
                F1,
                branch.gmin,
                branch.gmax,
                glo,
                ghi;
                h = h,
                segbuf = segbuf,
            )
            f2 = integrate_bin(
                F2,
                branch.gmin,
                branch.gmax,
                glo,
                ghi;
                h = h,
                segbuf = segbuf,
            )

            # find which bin to dump in
            t1 = 0.5 * (T1(glo) + T1(ghi)) - t0 + t_source_disc
            t2 = 0.5 * (T2(glo) + T2(ghi)) - t0 + t_source_disc
            i1 = searchsortedfirst(t_grid, t1)
            i2 = searchsortedfirst(t_grid, t2)

            imax = lastindex(t_grid)
            if i1 < imax
                flux[j, i1] += f1 * θ
            end
            if i2 < imax
                flux[j, i2] += f2 * θ
            end
        end
    end
end

export integrate_lineprofile, g✶_to_g, g_to_g✶
