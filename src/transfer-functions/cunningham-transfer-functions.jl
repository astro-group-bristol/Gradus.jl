function _make_interpolation(g, f)
    DataInterpolations.LinearInterpolation(f, g)
end

function _adjust_extrema!(g)
    g[1] = zero(eltype(g))
    g[end] = one(eltype(g))
end

function _make_sorted_with_adjustments!(g1, f1, g2, f2)
    I1 = sortperm(g1)
    I2 = sortperm(g2)

    # sort all in ascending g
    g1 .= g1[I1]
    g2 .= g2[I2]
    f1 .= f1[I1]
    f2 .= f2[I2]

    # at low inclination angles, the range of g is so tight that the edges are
    # very problematic due to the g✶ * (1 - g✶) terms, sending the branches to zero
    # so we use an offset to avoid these that is small enough to not impact results at
    # high inclination angles
    H = 1e-6
    J1 = @. (g1 < 1 - H) & (g1 > H)
    J2 = @. (g2 < 1 - H) & (g2 > H)
    g1 = g1[J1]
    g2 = g2[J2]
    f1 = f1[J1]
    f2 = f2[J2]

    _adjust_extrema!(g1)
    _adjust_extrema!(g2)
    g1, f1, g2, f2
end

function splitbranches(ctf::CunninghamTransferFunction{T}) where {T}
    # first things first we want to pull out the extremal
    # as these values are mutual to both branches, and cap off the
    # extrema. this way we avoid any accidental linear extrapolation
    # which might occur if N is small

    _, imin = findmin(ctf.g✶)
    _, imax = findmax(ctf.g✶)
    i1, i2 = imax > imin ? (imin, imax) : (imax, imin)

    if (i1 == i2)
        error("Resolved same min/max for rₑ = $(ctf.rₑ)")
    end

    # branch sizes, with duplicate extrema
    N1 = i2 - i1 + 1
    N2 = length(ctf.f) - N1 + 2
    # allocate branches
    branch1_f = zeros(T, N1)
    branch2_f = zeros(T, N2)
    branch1_g✶ = zeros(T, N1)
    branch2_g✶ = zeros(T, N2)

    for (i, j) in enumerate(i1:i2)
        branch1_f[i] = ctf.f[j]
        branch1_g✶[i] = ctf.g✶[j]
    end
    for (i, j) in enumerate(Iterators.flatten((1:i1, i2:length(ctf.f))))
        branch2_f[i] = ctf.f[j]
        branch2_g✶[i] = ctf.g✶[j]
    end

    # determine which is the upper branch
    # return:  (lower_g✶, lower_f, upper_g✶, upper_f)
    if branch1_f[2] > branch1_f[1]
        branch2_g✶, branch2_f, branch1_g✶, branch1_f
    else
        branch1_g✶, branch1_f, branch2_g✶, branch2_f
    end
end

function interpolate_transfer_function(ctf::CunninghamTransferFunction{T}) where {T}
    (lower_g✶, lower_f, upper_g✶, upper_f) = splitbranches(ctf)
    (lower_g✶, lower_f, upper_g✶, upper_f) =
        _make_sorted_with_adjustments!(lower_g✶, lower_f, upper_g✶, upper_f)
    lower_branch = _make_interpolation(lower_g✶, lower_f)
    upper_branch = _make_interpolation(upper_g✶, upper_f)
    InterpolatedCunninghamTransferFunction(
        upper_branch,
        lower_branch,
        ctf.gmin,
        ctf.gmax,
        ctf.rₑ,
    )
end

function _calculate_transfer_function(rₑ, g, g✶, J)
    if abs(g✶) > 1
        error("abs(g✶) > 1 for rₑ = $rₑ (g = $g, g✶ = $g✶, J = $J).")
    end
    @. (1 / (π * rₑ)) * g * √(g✶ * (1 - g✶)) * J
end

function cunningham_transfer_function(
    m::AbstractMetric{T},
    u,
    d,
    rₑ;
    max_time = 2 * u[2],
    redshift_pf = ConstPointFunctions.redshift(m, u),
    offset_max = 0.4rₑ + 10,
    zero_atol = 1e-7,
    θ_offset = 0.6,
    N = 80,
    tracer_kwargs...,
) where {T}
    function _workhorse(θ)
        r, gp = find_offset_for_radius(
            m,
            u,
            d,
            rₑ,
            θ;
            zero_atol = zero_atol,
            offset_max = offset_max,
            max_time = max_time,
            tracer_kwargs...,
        )
        if isnan(r)
            error(
                "Transfer function integration failed (rₑ=$rₑ, θ=$θ, offset_max = $offset_max).",
            )
        end
        α = r * cos(θ)
        β = r * sin(θ)
        # use underscores to avoid possible boxing
        _g = redshift_pf(m, gp, max_time)
        _J = jacobian_∂αβ_∂gr(
            m,
            u,
            d,
            α,
            β,
            max_time;
            redshift_pf = redshift_pf,
            tracer_kwargs...,
        )
        (_g, _J)
    end

    # sample so that the expected minima and maxima (0 and π)
    # are in the middle of the domain, so that we can find the minima
    # and maxima via interpolation
    # resample over domains of expected extrema to improve convergence
    M = (N ÷ 5)
    θs =
        Iterators.flatten((
            range(-π / 2, 3π / 2, N - 2 * M),
            range(-θ_offset, θ_offset, M),
            range(π - θ_offset, π + θ_offset, M),
        )) |> collect
    sort!(θs)

    # avoid coordinate singularity at π / 2
    @. θs += 1e-4

    Js = zeros(T, N)
    gs = zeros(T, N)
    @inbounds for i in eachindex(θs)
        θ = θs[i]
        g, J = _workhorse(θ)
        gs[i] = g
        Js[i] = J
    end

    # todo: maybe use a basic binary search optimization here to find
    # the extremal g using the ray tracer within N steps, instead of
    # relying on an interpolation to find it?
    # maybe even dispatch for different methods until one proves itself
    # superior

    # gmin, gmax = _search_extremal!(gs, Js, _workhorse, θs, Ndirect, Nextrema_solving)

    # interpolate the extrema
    _gmin, _gmax = try
        (_a, _b), _ = infer_extremal(gs, θs, 0, π; offset = θ_offset)
        _a, _b
    catch e
        if e isa Roots.ConvergenceFailed
            @warn ("Root finder failed to infer minima for rₑ = $rₑ. Using array extremal.")
            extrema(gs)
        else
            throw(e)
        end
    end
    gmin, gmax = _check_gmin_gmax(_gmin, _gmax, rₑ, gs)

    # convert from ∂g to ∂g✶
    @. Js = (gmax - gmin) * Js
    @inbounds for i in eachindex(Js)
        g✶ = g_to_g✶(gs[i], gmin, gmax)
        # Js is now storing f
        Js[i] = _calculate_transfer_function(rₑ, gs[i], g✶, Js[i])
        # gs is now storing g✶
        gs[i] = g✶
    end

    CunninghamTransferFunction(gs, Js, gmin, gmax, rₑ)
end

function _check_gmin_gmax(_gmin, _gmax, rₑ, gs)
    gmin = _gmin
    gmax = _gmax
    # sometimes the interpolation seems to fail if there are duplicate knots
    if isnan(gmin)
        @warn "gmin is NaN for rₑ = $rₑ (gmin = $gmin, extrema(gs) = $(extrema(gs)). Using extrema."
        gmin = minimum(gs)
    end
    if isnan(gmax)
        @warn "gmax is NaN for rₑ = $rₑ (gmax = $gmax, extrema(gs) = $(extrema(gs)). Using extrema."
        gmax = maximum(gs)
    end
    if gmin == gmax
        @warn (
            "gmin == gmax for rₑ = $rₑ (gmin = gmax = $gmin, extrema(gs) = $(extrema(gs)). Using extrema."
        )
        gmin, gmax = extrema(gs)
        if gmin == gmax
            error("Cannot use extrema")
        end
    end
    if gmin > minimum(gs)
        @warn (
            "Inferred minima > array minimum rₑ = $rₑ (gmin = $gmin, extrema(gs) = $(extrema(gs)). Using minimum."
        )
        gmin = minimum(gs)
    end
    if gmax < maximum(gs)
        @warn (
            "Inferred maximum < array maximum rₑ = $rₑ (gmax = $gmax, extrema(gs) = $(extrema(gs)). Using maximum."
        )
        gmax = maximum(gs)
    end
    return gmin, gmax
end

# find gmin and gmax via GoldenSection
# storing all attempts in gs and Js
function _search_extremal!(gs, Js, f, θs, offset, N)
    jmin = offset
    function _min_searcher(θ)
        @assert jmin <= N + offset
        g, Js[jmin] = f(θ)
        gs[jmin] = g
        jmin + 1
        g
    end

    jmax = offset + N
    function _max_searcher(θ)
        @assert jmax <= 2N + offset
        g, Js[jmax] = f(θ)
        gs[jmax] = g
        jmax + 1
        # invert 
        -g
    end

    _, imin = findmin(@view(gs[1:offset]))
    _, imax = findmax(@view(gs[1:offset]))

    # stride either side of our best guess so far
    res_min = Optim.optimize(
        _min_searcher,
        θs[imin-1],
        θs[imin+1],
        GoldenSection(),
        iterations = N,
    )
    res_max = Optim.optimize(
        _max_searcher,
        θs[imax-1],
        θs[imax+1],
        GoldenSection(),
        iterations = N,
    )
    # unpack result, remembering that maximum is inverted
    Optim.minimum(res_min), -Optim.minimum(res_max)
end

function infer_extremal(y, x, x0, x1; offset = 0.4)
    x0, y0 = interpolate_extremal(y, x, x0; offset = offset)
    x1, y1 = interpolate_extremal(y, x, x1; offset = offset)
    if y0 > y1
        return (y1, y0), (x1, x0)
    else
        return (y0, y1), (x0, x1)
    end
end

∂(f) = x -> ForwardDiff.derivative(f, x)
function interpolate_extremal(y, x, x0; offset = 0.4)
    interp = DataInterpolations.CubicSpline(y, x)
    x̄ = find_zero(∂(interp), (x0 - offset, x0 + offset), Roots.AlefeldPotraShi())
    x̄, interp(x̄)
end

function interpolated_transfer_branches(
    m::AbstractMetric{T},
    x,
    d,
    radii;
    verbose = false,
    kwargs...,
) where {T}
    progress_bar = init_progress_bar("Transfer functions:", length(radii), verbose)
    # IILF for calculating the interpolated branches
    𝔉 =
        rₑ -> begin
            # want to scale the initial position with radius
            # since redshift / jacobian values calculated at large impact parameters
            # seem to be inaccurate? either that or the root finder is up to something
            # but the problems seem to disappear by just keeping everything at low impact
            x_prob = SVector{4}(x[1], x[2], x[3], x[4])
            ctf = cunningham_transfer_function(
                m,
                x_prob,
                d,
                rₑ,
                ;
                chart = chart_for_metric(m, 10 * x_prob[2]),
                max_time = 10 * x_prob[2],
                kwargs...,
            )
            itp = interpolate_transfer_function(ctf)
            ProgressMeter.next!(progress_bar)
            itp
        end

    # calculate interpolated transfer functions for each emission radius
    _threaded_map(𝔉, radii)
end

export CunninghamTransferFunction,
    InterpolatedCunninghamTransferFunction,
    splitbranches,
    interpolate_transfer_function,
    cunningham_transfer_function
