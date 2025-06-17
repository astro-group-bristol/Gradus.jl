_promote_disc_for_transfer_functions(d::AbstractAccretionDisc) = (d, d)
function _promote_disc_for_transfer_functions(::ThinDisc{T}) where {T}
    plane = DatumPlane(zero(T))
    plane, plane
end

struct _TransferFunctionSetup{T,A}
    h::T
    θ_offset::T
    "Tolerance for root finding"
    zero_atol::T
    "Unstable radius with respect to when finer tolerances are needed for the Jacobian evaluation"
    unstable_radius::T
    α₀::T
    β₀::T
    N::Int
    N_extrema::Int
    # TODO: remove me
    root_solver::A
end

function _TransferFunctionSetup(
    m::AbstractMetric{T},
    d::AbstractAccretionGeometry;
    θ_offset = T(0.6),
    zero_atol = T(1e-7),
    N = 80,
    N_extrema = 17,
    α₀ = 0,
    β₀ = 0,
    h = T(1e-6),
    root_solver = nothing,
    kwargs...,
) where {T}
    # specialize the algorithm depending on whether we are calculating for a thin or thick disc
    _alg = if isnothing(root_solver)
        if d isa AbstractThickAccretionDisc
            RootsAlg()
        else
            NonLinearAlg()
        end
    else
        root_solver
    end
    setup = _TransferFunctionSetup{T,typeof(_alg)}(
        h,
        θ_offset,
        zero_atol,
        Gradus.isco(m) + 1,
        convert(T, α₀),
        convert(T, β₀),
        N,
        N_extrema,
        _alg,
    )
    kwargs, setup
end

_is_unstable(setup::_TransferFunctionSetup, d::AbstractAccretionDisc, r) =
    (d isa AbstractThickAccretionDisc) && (r < setup.unstable_radius)

_calculate_transfer_function(rₑ, g, g✶, J) = @. (1 / (π * rₑ)) * g * √(g✶ * (1 - g✶)) * J

function _adjust_extrema!(arr::AbstractArray{T}, min = zero(T), max = one(T)) where {T}
    arr[1] = min
    arr[end] = max
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
    _adjust_extrema!(t1, t_lo, t_hi)
    _adjust_extrema!(t2, t_lo, t_hi)
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
    TransferBranches{false}(
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
            zero_atol = setup.zero_atol,
            max_time = max_time,
            β₀ = setup.β₀,
            α₀ = setup.α₀,
            tracer_kwargs...,
        )
        if isnan(r)
            error("Transfer function integration failed (rₑ=$rₑ, θ=$θ).")
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
    disc, jacobian_disc = _promote_disc_for_transfer_functions(d)
    workhorse, jacobian_function, _ = _setup_workhorse_jacobian_with_kwargs(
        setup,
        m,
        x,
        disc,
        rₑ;
        jacobian_disc = jacobian_disc,
        kwargs...,
    )
    function _disc_workhorse(θ::T; whkwargs...)::NTuple{3,T} where {T}
        g, gp, r = workhorse(θ; whkwargs...)
        α, β = _rθ_to_αβ(r, θ; α₀ = setup.α₀, β₀ = setup.β₀)
        (g, jacobian_function(α, β), gp.x[1])
    end
end

function _rear_workhorse(
    setup::_TransferFunctionSetup,
    m::AbstractMetric,
    x,
    d::AbstractThickAccretionDisc,
    rₑ;
    max_time = 2 * x[2],
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
            jacobian_disc = d,
            kwargs...,
        )
    function _thick_workhorse(θ::T; whkwargs...)::NTuple{4,T} where {T}
        g, gp, r = datum_workhorse(θ; whkwargs...)
        r₊ = try
            r_thick, _ = _find_offset_for_radius(
                m,
                x,
                d,
                rₑ,
                θ;
                initial_r = r,
                zero_atol = setup.zero_atol,
                max_time = max_time,
                β₀ = setup.β₀,
                α₀ = setup.α₀,
                tracer_kwargs...,
                # don't echo warnings
                warn = false,
            )
            r_thick
        catch
            # if we fail, for whatever reason, to root solve on the thick discs,
            # we don't care, we just need a NaN value and then set that point to
            # "not visible"
            NaN
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
    data::_TransferDataAccumulator{T},
    setup::_TransferFunctionSetup,
    workhorse,
    θiterator,
    rₑ,
) where {T}
    θ_offset = setup.θ_offset
    for (i, θ) in enumerate(θiterator)
        insert_data!(data, i, θ, workhorse(θ))
    end

    # TODO: find where the largest difference in g is and refine

    N, gmin_candidate, gmax_candidate = _search_extremal!(data, workhorse, θ_offset)
    gmin, gmax = @views _check_gmin_gmax(gmin_candidate, gmax_candidate, rₑ, data.gs[1:N])
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
    solver_kwargs, setup = _TransferFunctionSetup(m, d; kwargs...)
    cunningham_transfer_function(setup, m, x, d, rₑ; solver_kwargs...)
end

function cunningham_transfer_function(
    setup::_TransferFunctionSetup,
    m::AbstractMetric,
    x,
    d::AbstractAccretionDisc,
    rₑ::T,
    ;
    chart = chart_for_metric(m, 2 * x[2]),
    max_time = 2 * x[2],
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
    gmin, gmax = _cunningham_transfer_function!(data, setup, workhorse, θiterator, rₑ)
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

    _, min_i = @views findmin(data.gs[1:i])
    _, max_i = @views findmax(data.gs[1:i])

    function _gmin_finder(θ)
        if i >= lastindex(data)
            error("i >= lastindex(data): $i >= $(lastindex(data))")
        end
        i += 1
        insert_data!(data, i, θ, workhorse(θ))
        data.gs[i]
    end
    function _gmax_finder(θ)
        -_gmin_finder(θ)
    end

    # stride either side of our best guess so far
    res_min = bracket_minimize(
        _gmin_finder,
        data.θs[min_i] - offset,
        data.θs[min_i] + offset,
        N = N,
    )
    res_max = bracket_minimize(
        _gmax_finder,
        data.θs[max_i] - offset,
        data.θs[max_i] + offset,
        N = N,
    )

    # unpack result, remembering that maximum is inverted
    i, res_min, -res_max
end

function bracket_minimize(f, low, high; N)
    y = one(low)
    a = low
    b = high
    for _ = 1:(N÷2)
        c = b - (b - a) * INVPHI
        d = a + (b - a) * INVPHI
        fc = f(c)
        fd = f(d)
        if isapprox(fc, fd, atol = 1e-7)
            break
        end

        if fc < fd
            y = min(y, fc)
            b = d
        else
            y = min(y, fd)
            a = c
        end
    end
    y
end

function interpolated_transfer_branches(
    m::AbstractMetric,
    x,
    d::AbstractAccretionDisc,
    radii;
    kwargs...,
)
    solver_kwargs, setup = _TransferFunctionSetup(m, d; kwargs...)
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

function transfer_function_grid(
    m::AbstractMetric,
    x,
    d::AbstractAccretionGeometry,
    radii;
    Ng = 20,
    kwargs...,
)
    itfs = interpolated_transfer_branches(m, x, d, radii; kwargs...)
    transfer_function_grid(itfs, Ng)
end

function transfer_function_grid(itfs::InterpolatingTransferBranches, Ng::Int)
    g✶_grid = collect(range(0, 1, Ng))
    r_grid = itfs.radii

    lower_f = [branch.lower_f(g) for g in g✶_grid, branch in itfs.branches]
    upper_f = [branch.upper_f(g) for g in g✶_grid, branch in itfs.branches]
    lower_t = [branch.lower_t(g) for g in g✶_grid, branch in itfs.branches]
    upper_t = [branch.upper_t(g) for g in g✶_grid, branch in itfs.branches]

    CunninghamTransferGrid(
        r_grid,
        g✶_grid,
        itfs.gmin,
        itfs.gmax,
        lower_f,
        upper_f,
        lower_t,
        upper_t,
    )
end

function make_transfer_function_table(
    M::Type{<:KerrMetric},
    d::AbstractAccretionDisc,
    a_range::AbstractVector,
    θ_range::AbstractVector;
    verbose = true,
    r_max = 500.0,
    n_radii = 150,
    kwargs...,
)
    function _mapper(a, θ)
        if verbose
            @info "a, θ = $((a, θ))"
        end
        m = M(1.0, a)
        x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)

        radii = Gradus.Grids._inverse_grid(Gradus.isco(m) + 1e-2, r_max, n_radii)
        Gradus.transfer_function_grid(m, x, d, radii; verbose = verbose, kwargs...)
    end
    grids = [_mapper(a, θ) for a in a_range, θ in θ_range]
    Gradus.CunninghamTransferTable((a_range, θ_range), grids)
end

"""
    transferfunctions(
        m::AbstractMetric,
        x,
        d::AbstractAccretionDisc;
        minrₑ = isco(m) + 1e-2,
        maxrₑ = 50,
        numrₑ = 100,
        radii = Grids._inverse_grid(minrₑ, maxrₑ, numrₑ),
        kwargs...,
    )

Pre-compute Cunningham transfer functions for the given model components.
Returns an [`InterpolatingTransferBranches`](@ref) instance.
"""
function transferfunctions(
    m::AbstractMetric,
    x,
    d::AbstractAccretionDisc;
    minrₑ = isco(m) + 1e-2,
    maxrₑ = 50,
    numrₑ = 100,
    radii = Grids._inverse_grid(minrₑ, maxrₑ, numrₑ),
    kwargs...,
)
    interpolated_transfer_branches(m, x, d, radii; kwargs...)
end

export CunninghamTransferData,
    TransferBranches,
    InterpolatingTransferBranches,
    splitbranches,
    interpolate_branches,
    cunningham_transfer_function,
    make_transfer_function_table,
    transferfunctions
