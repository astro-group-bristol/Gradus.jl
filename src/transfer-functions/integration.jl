_lineprofile_integrand(_, g) = g^2

struct IntegrationSetup{T,K,F,R,SegBufType}
    h::T
    time::K
    integrand::F
    pure_radial::R
    segbuf::SegBufType
end

IntegrationSetup(T::Type, integrand, pure_radial; time = nothing, h = 1e-8) =
    IntegrationSetup(h, time, integrand, pure_radial, alloc_segbuf(T))

function integrand(setup::IntegrationSetup, f, r, g, g✶)
    # the radial term and π is hoisted into the integration weight for that annulus
    ω = f * g / (√(g✶ * (1 - g✶)))
    ω * setup.integrand(r, g)
end

struct _IntegrationClosures{S<:IntegrationSetup,B<:TransferBranches}
    setup::S
    branch::B
end

function _time_bins(c::_IntegrationClosures, glo, ghi)
    g✶lo = clamp(g_to_g✶(glo, c.branch.gmin, c.branch.gmax), c.setup.h, 1 - c.setup.h)
    g✶hi = clamp(g_to_g✶(ghi, c.branch.gmin, c.branch.gmax), c.setup.h, 1 - c.setup.h)
    t_lo = (c.branch.lower_t(g✶lo) + c.branch.lower_t(g✶hi)) / 2
    t_hi = (c.branch.upper_t(g✶lo) + c.branch.upper_t(g✶hi)) / 2
    (t_lo, t_hi)
end

function _upper_branch(c::_IntegrationClosures)
    function _upper_branch_integrand(g)
        g✶ = g_to_g✶(g, c.branch.gmin, c.branch.gmax)
        f = c.branch.upper_f(g✶)
        integrand(c.setup, f, c.branch.rₑ, g, g✶)
    end
end

function _lower_branch(c::_IntegrationClosures)
    function _lower_branch_integrand(g)
        g✶ = g_to_g✶(g, c.branch.gmin, c.branch.gmax)
        f = c.branch.lower_f(g✶)
        integrand(c.setup, f, c.branch.rₑ, g, g✶)
    end
end

function _both_branches(c::_IntegrationClosures)
    function _both_branches_integrand(g)
        g✶ = g_to_g✶(g, c.branch.gmin, c.branch.gmax)
        f = c.branch.upper_f(g✶) + c.branch.lower_f(g✶)
        integrand(c.setup, f, c.branch.rₑ, g, g✶)
    end
end

g_to_g✶(g, gmin, gmax) = @. (g - gmin) / (gmax - gmin)
g✶_to_g(g✶, gmin, gmax) = @. (gmax - gmin) * g✶ + gmin

function integrate_edge(S, lim, lim_g✶, gmin, gmax, h)
    gh = g✶_to_g(lim_g✶, gmin, gmax)
    2 * S(gh) * abs(√gh - √lim) * (gmax - gmin) * √h
end

function integrate_bin(c::_IntegrationClosures, S, lo::T, hi) where {T}
    # alias for convenience
    gmin = c.branch.gmin
    gmax = c.branch.gmax
    h = c.setup.h

    glo = clamp(lo, gmin, gmax)
    ghi = clamp(hi, gmin, gmax)

    lum = zero(T)
    if glo == ghi
        return lum
    end

    g✶lo = g_to_g✶(lo, gmin, gmax)
    g✶hi = g_to_g✶(hi, gmin, gmax)

    if g✶lo < h
        if g✶hi > h
            lum += integrate_edge(S, glo, h, gmin, gmax, h)
            glo = g✶_to_g(h, gmin, gmax)
        else
            return integrate_edge(S, glo, g✶hi, gmin, gmax, h)
        end
    end

    if g✶hi > 1 - h
        if g✶lo < 1 - h
            lum += integrate_edge(S, ghi, 1 - h, gmin, gmax, h)
            ghi = g✶_to_g(1 - h, gmin, gmax)
        else
            return integrate_edge(S, ghi, g✶lo, gmin, gmax, h)
        end
    end

    res, _ = quadgk(S, glo, ghi; segbuf = c.setup.segbuf)
    lum += res
    return lum
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
        @. flux[i, :] = flux[i, :] / ḡ
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
    setup = IntegrationSetup(T, _lineprofile_integrand, ε; kwargs...)
    _integrate_transfer_problem!(flux, setup, itb, fine_rₑ_grid, g_grid)
    _normalize(flux, g_grid)
end

function integrate_lagtransfer(
    profile::AbstractDiscProfile,
    itb::InterpolatingTransferBranches{T},
    radii,
    g_grid,
    t_grid;
    Nr = 1000,
    t0 = 0,
    kwargs...,
) where {T}
    # pre-allocate output
    flux = zeros(T, (length(g_grid), length(t_grid)))
    fine_rₑ_grid = Grids._geometric_grid(extrema(radii)..., Nr) |> collect
    setup = IntegrationSetup(
        T,
        _lineprofile_integrand,
        profile.f.ε;
        time = r -> profile.t.t(r) - t0,
        kwargs...,
    )
    _integrate_transfer_problem!(flux, setup, itb, fine_rₑ_grid, g_grid, t_grid)
    _normalize(flux, g_grid)
end

function _integrate_transfer_problem!(
    output::AbstractVector,
    setup::IntegrationSetup{T,Nothing},
    itb::InterpolatingTransferBranches{T},
    radii_grid,
    g_grid,
) where {T}
    g_grid_view = @views g_grid[1:end-1]
    # build fine radial grid for trapezoidal integration
    @inbounds for (i, rₑ) in enumerate(radii_grid)
        closures = _IntegrationClosures(setup, itb(rₑ))
        Δrₑ = _trapezoidal_weight(radii_grid, rₑ, i)
        # integration weight for this annulus
        θ =
            Δrₑ * rₑ * setup.pure_radial(rₑ) * π /
            (closures.branch.gmax - closures.branch.gmin)
        @inbounds for j in eachindex(g_grid_view)
            glo = g_grid[j]
            ghi = g_grid[j+1]
            k = integrate_bin(closures, _both_branches(closures), glo, ghi)
            output[j] += k * θ
        end
    end
end

function _integrate_transfer_problem!(
    output::AbstractMatrix,
    setup::IntegrationSetup,
    itb::InterpolatingTransferBranches{T},
    radii_grid,
    g_grid,
    t_grid,
) where {T}
    g_grid_view = @views g_grid[1:end-1]
    # build fine radial grid for trapezoidal integration
    @inbounds for (i, rₑ) in enumerate(radii_grid)
        closures = _IntegrationClosures(setup, itb(rₑ))
        Δrₑ = _trapezoidal_weight(radii_grid, rₑ, i)
        θ =
            Δrₑ * rₑ * setup.pure_radial(rₑ) * π /
            (closures.branch.gmax - closures.branch.gmin)
        # time delay for this annuli
        t_source_disc = setup.time(rₑ)

        @inbounds for j in eachindex(g_grid_view)
            glo = g_grid[j]
            ghi = g_grid[j+1]
            k1 = integrate_bin(closures, _lower_branch(closures), glo, ghi)
            k2 = integrate_bin(closures, _upper_branch(closures), glo, ghi)

            # find which bin to dump in
            t_lo, t_hi = _time_bins(closures, glo, ghi)
            t1 = t_lo + t_source_disc
            t2 = t_hi + t_source_disc
            i1 = searchsortedfirst(t_grid, t1)
            i2 = searchsortedfirst(t_grid, t2)

            imax = lastindex(t_grid)
            if i1 <= imax
                output[j, i1] += k1 * θ
            end
            if i2 <= imax
                output[j, i2] += k2 * θ
            end
        end
    end
end

export integrate_lineprofile, g✶_to_g, g_to_g✶
