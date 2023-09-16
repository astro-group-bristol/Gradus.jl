struct _TransferFunctionSetup{T}
    h::T
    θ_offset::T
    "Unstable radius with respect to when finer tolerances are needed for the Jacobian evaluation"
    unstable_radius::T
    α₀::T
    β₀::T
    N::Int
    N_extrema::Int
end

function _TransferFunctionSetup(
    m::AbstractMetric{T};
    θ_offset = T(0.6),
    N = 80,
    N_extrema = 17,
    α₀ = 0,
    β₀ = 0,
    h = T(1e-6),
    kwargs...,
) where {T}
    setup = _TransferFunctionSetup{T}(
        h,
        θ_offset,
        Gradus.isco(m) + 1,
        convert(T, α₀),
        convert(T, β₀),
        N,
        N_extrema,
    )
    kwargs, setup
end

_is_unstable(setup::_TransferFunctionSetup, d::AbstractAccretionDisc, r) =
    (d isa AbstractThickAccretionDisc) && (r < setup.unstable_radius)

_calculate_transfer_function(rₑ, g, g✶, J) = @. (1 / (π * rₑ)) * g * √(g✶ * (1 - g✶)) * J

function _adjust_extrema!(g::AbstractArray{T}) where {T}
    g[1] = zero(T)
    g[end] = one(T)
end

function _make_sorted_with_adjustments!(g1, f1, t1, g2, f2, t2; h = 1e-6)
    I1 = sortperm(g1)
    I2 = sortperm(g2)

    Base.permute!(g1, I1)
    Base.permute!(f1, I1)
    Base.permute!(t1, I1)

    Base.permute!(g2, I2)
    Base.permute!(f2, I2)
    Base.permute!(t2, I2)

    # at low inclination angles, the range of g is so tight that the edges are
    # very problematic due to the g✶ * (1 - g✶) terms, sending the branches to zero
    # so we use an offset to avoid these that is small enough to not impact results at
    # high inclination angles
    J1 = @. (g1 < 1 - h) & (g1 > h)
    J2 = @. (g2 < 1 - h) & (g2 > h)
    g1 = g1[J1]
    g2 = g2[J2]
    f1 = f1[J1]
    f2 = f2[J2]

    # save endpoints in time, since these will be pretty accurate
    t_lo = (t1[1] + t2[1]) / 2
    t_hi = (t1[end] + t2[end]) / 2
    # mask
    t1 = t1[J1]
    t2 = t2[J2]
    # restore endpoints
    t1[1] = t_lo
    t1[end] = t_hi
    t2[1] = t_lo
    t2[end] = t_hi

    _adjust_extrema!(g1)
    _adjust_extrema!(g2)
    g1, f1, t1, g2, f2, t2
end

function splitbranches(ctf::CunninghamTransferData{T}) where {T}
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
    # todo: do this all at once with a matrix
    branch1_f = zeros(T, N1)
    branch2_f = zeros(T, N2)
    branch1_t = zeros(T, N1)
    branch2_t = zeros(T, N2)
    branch1_g✶ = zeros(T, N1)
    branch2_g✶ = zeros(T, N2)

    for (i, j) in enumerate(i1:i2)
        branch1_f[i] = ctf.f[j]
        branch1_g✶[i] = ctf.g✶[j]
        branch1_t[i] = ctf.t[j]
    end
    for (i, j) in enumerate(Iterators.flatten((1:i1, i2:length(ctf.f))))
        branch2_f[i] = ctf.f[j]
        branch2_g✶[i] = ctf.g✶[j]
        branch2_t[i] = ctf.t[j]
    end

    # determine which is the upper branch
    # return:  (lower_g✶, lower_f, upper_g✶, upper_f)
    if branch1_f[2] > branch1_f[1]
        branch2_g✶, branch2_f, branch2_t, branch1_g✶, branch1_f, branch1_t
    else
        branch1_g✶, branch1_f, branch1_t, branch2_g✶, branch2_f, branch2_t
    end
end

function interpolate_branches(ctf::CunninghamTransferData{T}; h = 1e-6) where {T}
    (lower_g✶, lower_f, lower_t, upper_g✶, upper_f, upper_t) = splitbranches(ctf)
    (lower_g✶, lower_f, lower_t, upper_g✶, upper_f, upper_t) =
        _make_sorted_with_adjustments!(
            lower_g✶,
            lower_f,
            lower_t,
            upper_g✶,
            upper_f,
            upper_t;
            h = h,
        )
    lower_branch = _make_interpolation(lower_g✶, lower_f)
    lower_time_branch = _make_interpolation(lower_g✶, lower_t)
    upper_branch = _make_interpolation(upper_g✶, upper_f)
    upper_time_branch = _make_interpolation(upper_g✶, upper_t)
    TransferBranches(
        upper_branch,
        lower_branch,
        upper_time_branch,
        lower_time_branch,
        ctf.gmin,
        ctf.gmax,
        ctf.rₑ,
    )
end

function _setup_workhorse_jacobian_with_kwargs(
    setup::_TransferFunctionSetup,
    m::AbstractMetric,
    x,
    d::AbstractAccretionDisc,
    rₑ;
    max_time = 2 * x[2],
    offset_max = 0.4rₑ + 10,
    zero_atol = 1e-8,
    redshift_pf = ConstPointFunctions.redshift(m, x),
    jacobian_disc = d,
    tracer_kwargs...,
)
    # underscores to avoid boxing variables
    function _jacobian_function(α_, β_)
        jacobian_∂αβ_∂gr(
            m,
            x,
            jacobian_disc,
            α_,
            β_,
            max_time;
            redshift_pf = redshift_pf,
            tracer_kwargs...,
        )
    end
    function _workhorse(θ)
        r, gp = find_offset_for_radius(
            m,
            x,
            d,
            rₑ,
            θ;
            zero_atol = zero_atol,
            offset_max = offset_max,
            max_time = max_time,
            β₀ = setup.β₀,
            α₀ = setup.α₀,
            tracer_kwargs...,
        )
        if isnan(r)
            error(
                "Transfer function integration failed (rₑ=$rₑ, θ=$θ, offset_max = $offset_max).",
            )
        end
        g = redshift_pf(m, gp, max_time)
        (g, gp, r)
    end
    (_workhorse, _jacobian_function, tracer_kwargs)
end
function _rear_workhorse(
    setup::_TransferFunctionSetup,
    m::AbstractMetric,
    x,
    d::AbstractAccretionDisc,
    rₑ;
    kwargs...,
)
    workhorse, jacobian_function, _ =
        _setup_workhorse_jacobian_with_kwargs(setup, m, x, d, rₑ; kwargs...)
    function _disc_workhorse(θ::T)::NTuple{3,T} where {T}
        g, gp, r = workhorse(θ)
        α, β = _rθ_to_αβ(r, θ; α₀ = setup.α₀, β₀ = setup.β₀)
        g, jacobian_function(α, β), gp.x[1]
    end
end

function _rear_workhorse(
    setup::_TransferFunctionSetup,
    m::AbstractMetric,
    x,
    d::AbstractThickAccretionDisc,
    rₑ;
    max_time = 2 * x[2],
    zero_atol = 1e-8,
    offset_max = 0.4rₑ + 10,
    kwargs...,
)
    plane = datumplane(d, rₑ)
    datum_workhorse, jacobian_function, tracer_kwargs =
        _setup_workhorse_jacobian_with_kwargs(
            setup,
            m,
            x,
            plane,
            rₑ;
            max_time = max_time,
            zero_atol = zero_atol,
            offset_max = offset_max,
            jacobian_disc = d,
            kwargs...,
        )
    function _thick_workhorse(θ::T)::NTuple{4,T} where {T}
        g, gp, r = datum_workhorse(θ)
        r₊, _ = try
            _find_offset_for_radius(
                m,
                x,
                d,
                rₑ,
                θ;
                initial_r = r,
                zero_atol = zero_atol,
                offset_max = offset_max,
                max_time = max_time,
                β₀ = setup.β₀,
                α₀ = setup.α₀,
                tracer_kwargs...,
                # don't echo warnings
                warn = false,
            )
        catch err
            if err isa Roots.ConvergenceFailed
                # return gp for type stability
                NaN, gp
            else
                throw(err)
            end
        end
        is_visible, J = if !isnan(r₊) && isapprox(r, r₊, atol = 1e-3)
            # trace jacobian on updated impact parameters
            α, β = _rθ_to_αβ(r₊, θ; α₀ = setup.α₀, β₀ = setup.β₀)
            _J = jacobian_function(α, β)
            isfinite(_J), _J
        else
            false, NaN
        end
        (g, J, gp.x[1], is_visible)
    end
end

function _cunningham_transfer_function!(
    data::_TransferDataAccumulator,
    workhorse,
    θiterator,
    θ_offset,
    rₑ,
)
    for (i, θ) in enumerate(θiterator)
        θ_corrected = θ + 1e-4
        insert_data!(data, i, θ_corrected, workhorse(θ))
    end

    gmin_candidate, gmax_candidate = _search_extremal!(data, workhorse, θ_offset)
    gmin, gmax = _check_gmin_gmax(gmin_candidate, gmax_candidate, rₑ, data.gs)
    # we might not have used all of the memory we allocated so let's clean up
    remove_unused_elements!(data)
    # sort everything by angle
    sort!(data)

    # convert from ∂g to ∂g✶
    @. data.Js = (gmax - gmin) * data.Js

    @inbounds for i in eachindex(data)
        g✶ = g_to_g✶(data.gs[i], gmin, gmax)
        # Js is now storing f
        data.Js[i] = _calculate_transfer_function(rₑ, data.gs[i], g✶, data.Js[i])
        # gs is now storing g✶
        data.gs[i] = g✶
    end

    gmin, gmax
end

function cunningham_transfer_function(
    m::AbstractMetric{Q},
    x,
    d::AbstractAccretionDisc,
    rₑ::T;
    kwargs...,
) where {Q,T}
    solver_kwargs, setup = _TransferFunctionSetup(m; kwargs...)
    cunningham_transfer_function(setup, m, x, d, rₑ; solver_kwargs...)
end

function cunningham_transfer_function(
    setup::_TransferFunctionSetup,
    m::AbstractMetric,
    x,
    d::AbstractAccretionDisc,
    rₑ::T,
    ;
    chart = chart_for_metric(m, 10 * x[2]),
    max_time = 10 * x[2],
    solver_kwargs...,
) where {T}
    M = setup.N + 2 * setup.N_extrema
    K = setup.N ÷ 5
    data = _TransferDataAccumulator(T, M, setup.N)
    # sample so that the expected minima and maxima (0 and π)
    θiterator = Iterators.flatten((
        range(0 - 2setup.θ_offset, 0 + 2setup.θ_offset, K),
        range(-π / 2, 3π / 2, setup.N - 2 * K),
        range(π - 2setup.θ_offset, π + 2setup.θ_offset, K),
    ))
    workhorse = _rear_workhorse(
        setup,
        m,
        x,
        d,
        rₑ;
        max_time = max_time,
        chart = chart,
        solver_kwargs...,
    )
    gmin, gmax =
        _cunningham_transfer_function!(data, workhorse, θiterator, setup.θ_offset, rₑ)
    CunninghamTransferData(
        data.data[2, :],
        data.data[3, :],
        data.data[4, :],
        gmin,
        gmax,
        rₑ,
    )
end

# find gmin and gmax via GoldenSection
# storing all attempts in gs and Js
function _search_extremal!(data::_TransferDataAccumulator, workhorse, offset)
    # need to specify the type to avoid boxing
    i::Int = data.cutoff
    N = (lastindex(data) - data.cutoff) ÷ 2 - 1

    function _gmin_finder(θ)
        if i >= lastindex(data)
            error("i >= lastindex(data): $i >= $(lastindex(data))")
        end
        i += 1
        if abs(θ) < 1e-4 || abs(abs(θ) - π) < 1e-4
            θ += 1e-4
        end
        insert_data!(data, i, θ, workhorse(θ))
        data.gs[i]
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
        iterations = N,
    )
    res_max = Optim.optimize(
        _gmax_finder,
        π - offset,
        π + offset,
        GoldenSection(),
        iterations = N,
    )

    # unpack result, remembering that maximum is inverted
    Optim.minimum(res_min), -Optim.minimum(res_max)
end

function interpolated_transfer_branches(
    m::AbstractMetric,
    x,
    d::AbstractAccretionDisc,
    radii;
    kwargs...,
)
    solver_kwargs, setup = _TransferFunctionSetup(m; kwargs...)
    interpolated_transfer_branches(setup, m, x, d, radii; solver_kwargs...)
end

function interpolated_transfer_branches(
    setup::_TransferFunctionSetup,
    m::AbstractMetric{T},
    x,
    d,
    radii;
    verbose = false,
    solver_kwargs...,
) where {T}
    progress_bar = init_progress_bar("Transfer functions:", length(radii), verbose)
    # IILF for calculating the interpolated branches
    𝔉 =
        rₑ -> begin
            ctf = cunningham_transfer_function(
                setup,
                m,
                x,
                d,
                rₑ,
                ;
                verbose = verbose,
                solver_kwargs...,
            )
            itp = interpolate_branches(ctf; h = setup.h)
            ProgressMeter.next!(progress_bar)
            itp
        end
    # calculate interpolated transfer functions for each emission radius
    InterpolatingTransferBranches(_threaded_map(𝔉, radii))
end

export CunninghamTransferData,
    TransferBranches,
    InterpolatingTransferBranches,
    splitbranches,
    interpolate_branches,
    cunningham_transfer_function
