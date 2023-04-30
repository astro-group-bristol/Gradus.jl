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
    m::AbstractMetric,
    u,
    d,
    rₑ::T;
    max_time = 2 * u[2],
    redshift_pf = ConstPointFunctions.redshift(m, u),
    offset_max = 0.4rₑ + 10,
    zero_atol = 1e-7,
    θ_offset = 0.6,
    N = 80,
    N_extrema = 17,
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

    M = N + 2 * N_extrema
    K = N ÷ 5

    θs = zeros(T, M)
    Js = zeros(T, M)
    gs = zeros(T, M)
    # sample so that the expected minima and maxima (0 and π)
    itt = Iterators.flatten((
        range(0 - 2θ_offset, 0 + 2θ_offset, K),
        range(-π / 2, 3π / 2, N - 2 * K),
        range(π - 2θ_offset, π + 2θ_offset, K),
    ))
    @inbounds for (i, _θ) in enumerate(itt)
        # avoid coordinate singularity at π / 2
        θ = _θ + 1e-4
        g, J = _workhorse(θ)
        gs[i] = g
        Js[i] = J
        θs[i] = θ
    end

    _gmin, _gmax = @views _search_extremal!(
        θs[N+1:end],
        gs[N+1:end],
        Js[N+1:end],
        _workhorse,
        θ_offset,
    )

    # remove unused 
    B = findall(==(0), θs)
    deleteat!(θs, B)
    deleteat!(gs, B)
    deleteat!(Js, B)

    gmin, gmax = _check_gmin_gmax(_gmin, _gmax, rₑ, gs)

    I = sortperm(θs)
    @. θs = θs[I]
    @. Js = Js[I]
    @. gs = gs[I]
    
    # convert from ∂g to ∂g✶
    @. Js = (gmax - gmin) * Js

    @inbounds for i in eachindex(gs)
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
function _search_extremal!(θs, gs, Js, f, offset)
    N = length(θs)

    i::Int = 0
    function _gmin_finder(θ)
        if i > N
            error("i > N: $i > $N")
        end
        i += 1
        if abs(θ) < 1e-4 || abs(abs(θ) - π) < 1e-4
            θ += 1e-4
        end
        g, J = f(θ)
        θs[i] = θ
        Js[i] = J
        gs[i] = g
        g
    end
    function _gmax_finder(θ)
        -_gmin_finder(θ)
    end

    # stride either side of our best guess so far
    res_min = Optim.optimize(
        _gmin_finder,
        0 - offset,
        0 + offset,
        GoldenSection(),
        iterations = N ÷ 2 - 1,
    )
    res_max = Optim.optimize(
        _gmax_finder,
        π - offset,
        π + offset,
        GoldenSection(),
        iterations = N ÷ 2 - 1,
    )

    # unpack result, remembering that maximum is inverted
    Optim.minimum(res_min), -Optim.minimum(res_max)
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
