_lineprofile_integrand(_, g) = g^2

function quadrature_integrate(f, a::T, b::T; rule = gauss(5)) where {T}
    total = zero(T)
    for (xᵢ, wᵢ) in zip(rule...)
        u = (xᵢ + 1) * ((b - a) / 2) + a
        total += wᵢ * f(u)
    end
    total * (b - a) / 2
end

struct IntegrationSetup{T,K,F,R}
    h::T
    time::K
    integrand::F
    pure_radial::R
    quadrature_rule::Tuple{Vector{T},Vector{T}}
    index_cache::Vector{Int}
    g_grid_upscale::Int
    n_radii::Int
end

function integration_setup(prof, transfer_functions; kwargs...)
    radial = if prof isa AbstractDiscProfile
        r -> emissivity_at(prof, r)
    else
        prof
    end
    IntegrationSetup(_lineprofile_integrand, radial; kwargs...)
end

function IntegrationSetup(
    integrand,
    pure_radial;
    time = nothing,
    h = 1e-8,
    g_grid_upscale = 1,
    n_radii = 1000,
    quadrature_points = 23,
)
    IntegrationSetup(
        h,
        time,
        integrand,
        pure_radial,
        gauss(quadrature_points),
        ones(Int, 4),
        g_grid_upscale,
        n_radii,
    )
end

function _g_fine_grid_iterate(upscale, lo, hi)
    Δ = (hi - lo) / upscale
    range(lo, hi - Δ, upscale), Δ
end

function integrand(setup::IntegrationSetup, f, r, g, g✶)
    # the radial term and π is hoisted into the integration weight for that annulus
    ω = f * g / (√(g✶ * (1 - g✶)))
    ω * setup.integrand(r, g)
end

function _time_interpolate(t0, f1, f2, g✶, h)
    ω, t1, t2 = if g✶ < h
        i = g✶ / h
        i, f1(h), f2(h)
    elseif g✶ > 1 - h
        i = 1 - (1 - g✶) / h
        i, f1(1 - h), f2(1 - h)
    else
        return t0
    end
    t1 * ω + (1 - ω) * t2
end

function _time_g✶(setup::IntegrationSetup, branch::TransferBranches, g✶)
    t_lower = branch.lower_t(g✶)
    t_upper = branch.upper_t(g✶)
    t1 = _time_interpolate(t_lower, branch.lower_t, branch.upper_t, g✶, setup.h)
    t2 = _time_interpolate(t_upper, branch.upper_t, branch.lower_t, g✶, setup.h)
    t1, t2
end

function _time_bins(setup::IntegrationSetup, branch::TransferBranches, glo, ghi)
    g✶lo = clamp(g_to_g✶(glo, branch.gmin, branch.gmax), 0, 1)
    g✶hi = clamp(g_to_g✶(ghi, branch.gmin, branch.gmax), 0, 1)

    tl1, tu1 = _time_g✶(setup, branch, g✶lo)
    tl2, tu2 = _time_g✶(setup, branch, g✶hi)
    (tl1, tl2), (tu1, tu2)
end

function _upper_branch(setup::IntegrationSetup, branch::TransferBranches)
    function _upper_branch_integrand(g)
        g✶ = g_to_g✶(g, branch.gmin, branch.gmax)
        f = branch.upper_f(g✶)
        integrand(setup, _zero_if_nan(f), branch.rₑ, g, g✶)
    end
end

function _lower_branch(setup::IntegrationSetup, branch::TransferBranches)
    function _lower_branch_integrand(g)
        g✶ = g_to_g✶(g, branch.gmin, branch.gmax)
        f = branch.lower_f(g✶)
        integrand(setup, _zero_if_nan(f), branch.rₑ, g, g✶)
    end
end

# find the index of the first `x` where `x[i] >= t`
function _search_sorted_chronological(x, t, last_i)
    for i = last_i:lastindex(x)
        if x[i] >= t
            return i
        end
    end
    return lastindex(x)
end

function _both_branches(
    setup::IntegrationSetup,
    branch::TransferBranches{SameDomain},
) where {SameDomain}
    function _both_branches_integrand(g)
        g✶ = g_to_g✶(g, branch.gmin, branch.gmax)

        fl, fu = if SameDomain
            x = branch.upper_f.f1.t

            i1 = _search_sorted_chronological(x, g✶, setup.index_cache[1])
            setup.index_cache[1] = i1
            weight = (g✶ - x[i1-1]) / (x[i1] - x[i1-1])

            _fu = branch.upper_f(i1 - 1, weight)
            _fl = branch.lower_f(i1 - 1, weight)
            (_fl, _fu)
        else
            _fu = branch.upper_f(g✶)
            _fl = branch.lower_f(g✶)
            (_fl, _fu)
        end

        f = _zero_if_nan(fu) + _zero_if_nan(fl)
        integrand(setup, f, branch.rₑ, g, g✶)
    end
end

g_to_g✶(g, gmin, gmax) = @. (g - gmin) / (gmax - gmin)
g✶_to_g(g✶, gmin, gmax) = @. (gmax - gmin) * g✶ + gmin

function integrate_edge(S, lim, lim_g✶, gmin, gmax)
    gh = g✶_to_g(lim_g✶, gmin, gmax)
    2 * S(gh) * abs(√gh - √lim)
end

# a very important `@inline`
@inline function integrate_bin(setup::IntegrationSetup, S, lo::T, hi, gmin, gmax) where {T}
    h = setup.h
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

    res = quadrature_integrate(S, glo, ghi; rule = setup.quadrature_rule)
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
    else
        x - X[i-1]
    end
end

function integrate_lineprofile(
    prof,
    transfer_functions,
    g_grid;
    rmin = inner_radius(transfer_functions),
    rmax = outer_radius(transfer_functions),
    kwargs...,
)
    setup = integration_setup(prof, transfer_functions; kwargs...)
    output = zeros(eltype(g_grid), length(g_grid))
    integrate_lineprofile!(
        output,
        setup,
        transfer_functions,
        g_grid;
        rmin = rmin,
        rmax = rmax,
    )
    output
end

function integrate_lagtransfer(
    prof,
    transfer_functions,
    g_grid,
    t_grid;
    rmin = inner_radius(transfer_functions),
    rmax = outer_radius(transfer_functions),
    t0 = 0,
    kwargs...,
)
    setup = integration_setup(
        prof,
        transfer_functions;
        time = r -> coordtime_at(prof, r) - t0,
        kwargs...,
    )
    output = zeros(eltype(g_grid), (length(g_grid), length(t_grid)))
    integrate_lagtransfer!(
        output,
        setup,
        transfer_functions,
        g_grid,
        t_grid;
        rmin = rmin,
        rmax = rmax,
    )
    output
end

function integrate_lineprofile!(
    output::AbstractVector,
    setup::IntegrationSetup,
    transfer_functions,
    grid::AbstractVector;
    rmin = inner_radius(transfer_functions),
    rmax = outer_radius(transfer_functions),
)
    # zero the output
    output .= zero(eltype(output))
    _integrate_transfer_problem!(output, setup, transfer_functions, (rmin, rmax), grid)
    _normalize!(output, grid)
end

function integrate_lagtransfer!(
    output::AbstractMatrix,
    setup::IntegrationSetup,
    transfer_functions,
    grid::AbstractVector,
    t_grid::AbstractVector;
    rmin = inner_radius(transfer_functions),
    rmax = outer_radius(transfer_functions),
)
    # zero the output
    output .= zero(eltype(output))
    _integrate_transfer_problem!(
        output,
        setup,
        transfer_functions,
        (rmin, rmax),
        grid,
        t_grid,
    )
    _normalize!(output, grid)
end

function _integrate_transfer_problem!(
    output::AbstractVector,
    setup::IntegrationSetup{T,Nothing},
    transfer_function_radial_interpolation,
    r_limits,
    g_grid,
) where {T}
    g_grid_view = @views g_grid[1:end-1]

    r_itterator = Grids._inverse_grid(r_limits..., setup.n_radii)
    r2 = first(iterate(r_itterator, 2))
    # prime the first r_prev so that the bin width is r2 - r1
    r_prev = r_limits[1] - (r2 - r_limits[1])

    for rₑ in r_itterator
        branch = transfer_function_radial_interpolation(rₑ)
        S = _both_branches(setup, branch)

        Δrₑ = rₑ - r_prev
        # integration weight for this annulus
        θ = Δrₑ * rₑ * setup.pure_radial(rₑ) * π / (branch.gmax - branch.gmin)

        @inbounds for j in eachindex(g_grid_view)
            glo = g_grid[j]
            ghi = g_grid[j+1]
            flux_contrib = zero(eltype(output))
            Δg = (ghi - glo) / setup.g_grid_upscale

            for i = 1:setup.g_grid_upscale
                g_fine_lo = glo + (i - 1) * Δg
                g_fine_hi = g_fine_lo + Δg
                k = integrate_bin(setup, S, g_fine_lo, g_fine_hi, branch.gmin, branch.gmax)
                flux_contrib += k
            end
            output[j] += flux_contrib * θ
        end

        r_prev = rₑ
    end
end

function _integrate_transfer_problem!(
    output::AbstractMatrix,
    setup::IntegrationSetup{T},
    transfer_function_radial_interpolation,
    r_limits,
    g_grid,
    t_grid,
) where {T}
    g_grid_view = @views g_grid[1:end-1]

    r_itterator = Grids._geometric_grid(r_limits..., setup.n_radii)
    r2 = first(iterate(r_itterator, 2))
    # prime the first r_prev so that the bin width is r2 - r1
    r_prev = r_limits[1] - (r2 - r_limits[1])

    for rₑ in r_itterator
        branch = transfer_function_radial_interpolation(rₑ)
        S_lower = _lower_branch(setup, branch)
        S_upper = _upper_branch(setup, branch)

        Δrₑ = rₑ - r_prev
        # integration weight for this annulus
        θ = Δrₑ * rₑ * setup.pure_radial(rₑ) * π / (branch.gmax - branch.gmin)

        # time delay for this annuli
        t_source_disc = setup.time(rₑ)

        @inbounds for j in eachindex(g_grid_view)
            glo = clamp(g_grid[j], branch.gmin, branch.gmax)
            ghi = clamp(g_grid[j+1], branch.gmin, branch.gmax)
            # skip if bin not relevant
            if glo == ghi
                continue
            end

            Δg = (ghi - glo) / setup.g_grid_upscale

            for i = 1:setup.g_grid_upscale
                g_fine_lo = glo + (i - 1) * Δg
                g_fine_hi = g_fine_lo + Δg

                k1 = integrate_bin(
                    setup,
                    S_lower,
                    g_fine_lo,
                    g_fine_hi,
                    branch.gmin,
                    branch.gmax,
                )
                k2 = integrate_bin(
                    setup,
                    S_upper,
                    g_fine_lo,
                    g_fine_hi,
                    branch.gmin,
                    branch.gmax,
                )

                # find which bin to dump in
                (tl1, tl2), (tu1, tu2) = _time_bins(setup, branch, g_fine_lo, g_fine_hi)
                i1 = searchsortedfirst(t_grid, (tl1 + tl2) / 2 + t_source_disc)
                i2 = searchsortedfirst(t_grid, (tu1 + tu2) / 2 + t_source_disc)

                imax = lastindex(t_grid)
                if i1 <= imax
                    output[j, i1] += k1 * θ
                end
                if i2 <= imax
                    output[j, i2] += k2 * θ
                end
            end
        end

        r_prev = rₑ
    end
end

export integrate_lineprofile, g✶_to_g, g_to_g✶, integration_setup, IntegrationSetup
