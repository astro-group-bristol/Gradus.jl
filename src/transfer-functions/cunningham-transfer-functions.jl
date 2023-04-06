function _make_sorted_interpolation(g, f)
    I = sortperm(g)
    _g = @inbounds g[I]
    _f = @inbounds f[I]
    # shift the extrema
    # todo: interpolate f -> g to find better estimate of what f should be at extrema
    _g[1] = zero(eltype(_g))
    _g[end] = one(eltype(_g))
    DataInterpolations.LinearInterpolation(_f, _g)
end

function splitbranches(ctf::CunninghamTransferFunction{T}) where {T}
    # first things first we want to pull out the extremal
    # as these values are mutual to both branches, and cap off the
    # extrema. this way we avoid any accidental linear extrapolation
    # which might occur if N is small

    _, imin = findmin(ctf.g✶)
    _, imax = findmax(ctf.g✶)
    i1, i2 = imax > imin ? (imin, imax) : (imax, imin)

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
    lower_branch = _make_sorted_interpolation(lower_g✶, lower_f)
    upper_branch = _make_sorted_interpolation(upper_g✶, upper_f)
    InterpolatedCunninghamTransferFunction(
        upper_branch,
        lower_branch,
        ctf.gmin,
        ctf.gmax,
        ctf.rₑ,
    )
end

_calculate_transfer_function(rₑ, g, g✶, J) = @. (1 / (π * rₑ)) * g * √(g✶ * (1 - g✶)) * J

function cunningham_transfer_function(
    m::AbstractMetric{T},
    u,
    d,
    rₑ;
    max_time = 2e3,
    diff_order = 4,
    redshift_pf = ConstPointFunctions.redshift(m, u),
    offset_max = 20.0,
    zero_atol = 1e-7,
    N = 80,
    tracer_kwargs...,
) where {T}
    # add 2 for extremal g✶ we calculate
    # Nextrema_solving = fld(N, 4)
    # Ndirect = N - 2 * (Nextrema_solving)
    Js = zeros(T, N)
    gs = zeros(T, N)

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
            error("Transfer function integration failed (rₑ=$rₑ, θ=$θ).")
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
            diff_order = diff_order,
            redshift_pf = redshift_pf,
            tracer_kwargs...,
        )
        (_g, _J)
    end

    # sample so that the expected minima and maxima (π and 2π)
    # are in the middle of the domain, so that we can find the minima
    # and maxima via interpolation
    θs = range(π / 2, 2π + π / 2, N) |> collect
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
    (gmin, gmax), _ = infer_extremal(gs, θs, π, 2π)

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

# currently unused: find gmin and gmax via GoldenSection
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

function infer_extremal(y, x, x0, x1)
    x0, y0 = interpolate_extremal(y, x, x0)
    x1, y1 = interpolate_extremal(y, x, x1)
    if y0 > y1
        return (y1, y0), (x1, x0)
    else
        return (y0, y1), (x0, x1)
    end
end

∂(f) = x -> ForwardDiff.derivative(f, x)
function interpolate_extremal(y, x, x0)
    interp = DataInterpolations.CubicSpline(y, x)
    x̄ = find_zero(∂(interp), x0)
    x̄, interp(x̄)
end

function interpolated_transfer_branches(
    m::AbstractMetric{T},
    u,
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
            u_prob = SVector{4}(u[1], 1000 + 100rₑ, u[3], u[4])
            ctf = cunningham_transfer_function(
                m,
                u_prob,
                d,
                rₑ,
                ;
                offset_max = 3rₑ + 20.0,
                chart = chart_for_metric(m, 10 * u_prob[2]),
                max_time = 10 * u_prob[2],
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
