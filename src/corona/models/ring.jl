const INVPHI = (sqrt(5) - 1) / 2

"""
A thread local cache for doing the calculations related to tracing the arms of
a point-approximation of a ring corona.

One is created, which is then copied once for each thread.
"""
struct _RingCoronaCache{T,S<:EmissivityProfileSetup{true},A}
    setup::S
    "Position of a point in the ring."
    x::SVector{4,T}
    "Velocity of the point in the ring."
    v::SVector{4,T}
    "(Local) poloidal angles that have been calculated."
    angles::Vector{T}
    "Geodesic points calculated"
    gps::Vector{GeodesicPoint{T,A}}
    "An optimising numeric buffer that is the same length of `gps` for storing
    partial results."
    _buffer::Vector{T}
    "Used to track which angles and gps we have written so far. `_index == 0` means the array is empty."
    _index::Ref{Int}
    "Number of geodesics for all the arms"
    N::Int
    "Maximal number of iterations that the extremiser should make"
    extrema_iter::Int
end

function Base.copy(cache::_RingCoronaCache)
    _RingCoronaCache(
        cache.setup,
        cache.x,
        cache.v,
        deepcopy(cache.angles),
        deepcopy(cache.gps),
        deepcopy(cache._buffer),
        Ref(cache._index[]),
        cache.N,
        cache.extrema_iter,
    )
end

function _add_to_cache!(
    cache::_RingCoronaCache{T,S,A},
    ang::T,
    gp::GeodesicPoint{T,A},
) where {T,S,A}
    cache._index[] += 1
    i::Int = cache._index[]

    if i > lastindex(cache.angles)
        # TODO: some kind of warning?
        @warn "Not enough pre-allocated space, ignoring values!"
        return
    end

    cache.angles[i] = ang
    cache.gps[i] = gp
end

function _RingCoronaCache(
    setup::EmissivityProfileSetup,
    m::AbstractMetric,
    model::RingCorona;
    h = 1e-7,
    extrema_iter = 80,
)
    x, v = sample_position_velocity(m, model)

    # TODO: infer this
    GPType = GeodesicPoint{eltype(x),Nothing}

    N = (setup.n_samples - 2 * extrema_iter)
    if (N <= 0)
        error(
            "Too few samples for selected method. Pass `n_samples` of more than $(2 * extrema_iter)",
        )
    end
    if !iseven(N)
        error("`n_samples` must be an even number for this method (given $N)")
    end

    angles = zeros(eltype(x), setup.n_samples)
    _buffer = zeros(eltype(x), setup.n_samples)
    gps = Vector{GPType}(undef, setup.n_samples)

    num_per_arm = div(N, 2)

    # setup the angles we are going to sample evenly over the sky
    left_arm_angles = range(h, π - h, num_per_arm)
    for (i, a) in enumerate(left_arm_angles)
        angles[i] = a
        angles[i+num_per_arm] = mod2pi(a + π)
    end
    sort!(@views(angles[1:N]))

    _RingCoronaCache(setup, x, v, angles, gps, _buffer, Ref(0), N, extrema_iter)
end

"""
Returns a function that takes position, velocity, local angle, and the slice
angle and returns the initial velocity vector for the geodesics.
"""
function rotatorfunctor(m::AbstractMetric{T}, x::SVector, v::SVector, β) where {T}
    angle_of_axis = x[3]
    k = Gradus._cart_local_direction(angle_of_axis, zero(T))

    function _rotator(θ)
        q = _cart_local_direction(θ + angle_of_axis, zero(T))
        b = rodrigues_rotate(k, q, β)
        # convert back to spherical coordinates
        ph = atan(b[2], b[1])
        th = atan(sqrt(b[2]^2 + b[1]^2), b[3])

        sky_angles_to_velocity(m, x, v, th, ph)
    end
end

function determine_bracket(
    N::Int,
    angles::Vector,
    gps::Vector{<:GeodesicPoint{T}},
    comparator::Function,
) where {T}
    a = first(angles)
    b = first(angles)
    estimate::Union{Nothing,T} = nothing
    for i = 1:N
        gp = gps[i]

        if gp.status == StatusCodes.IntersectedWithGeometry
            rho = _equatorial_project(gp.x)

            if isnothing(estimate) || comparator(rho, estimate)
                estimate = rho
                a = b = angles[i]

                if i < N && gps[i+1].status != StatusCodes.IntersectedWithGeometry
                    # this will likely be a better upper limit
                    b = angles[i+1]
                end

                if i > 1 && gps[i-1].status != StatusCodes.IntersectedWithGeometry
                    # this will likely be a better lower limit
                    a = angles[i-1]
                end
            end
        end
    end

    best_estimate::T = estimate
    delta = π / min(N, 32)
    a - 2 * delta, b + 2 * delta, best_estimate
end

"""
Given an already traced set of angles, use the best estimate of where the maxima
as and then bisect to try and find a better maxima.
"""
function _golden_bracket!(
    cache::_RingCoronaCache{T},
    tracer::Function,
    target::Val{Target};
    kwargs...,
) where {T,Target}
    comparator = if Target == :minima
        (a, b) -> a < b
    elseif Target == :maxima
        (a, b) -> a > b
    else
        error("Unknown target: $Target")
    end

    a, b, best_estimate = determine_bracket(cache.N, cache.angles, cache.gps, comparator)
    a_init, b_init = a, b

    function _objective(θ)
        gp = tracer(θ)
        if gp.status != StatusCodes.IntersectedWithGeometry
            # just return the best estimate, as this will (hopefully) be worse
            # than whatever we've calculated so far
            gp, best_estimate
        else
            gp, _equatorial_project(gp.x)
        end
    end

    # take a bracketing interval
    c_value = 0.0

    n = 0
    iters = 0
    too_many_iters = false
    while (n < cache.extrema_iter)
        c = b - (b - a) * INVPHI
        d = a + (b - a) * INVPHI

        c_gp, c_value = _objective(c)
        d_gp, d_value = _objective(d)

        if comparator(c_value, d_value)
            if comparator(c_value, best_estimate)
                n += 1
                _add_to_cache!(cache, c, c_gp)
            elseif too_many_iters
                n += 1
            end
            b = d
        else
            if comparator(d_value, best_estimate)
                n += 1
                _add_to_cache!(cache, d, d_gp)
            elseif too_many_iters
                n += 1
            end
            a = c
        end

        iters += 1
        if !too_many_iters && (iters > cache.extrema_iter)
            @warn "Too many iterations solving for $Target" maxlog = 1
            too_many_iters = true
        end
    end
end

"""
    function _ring_arm_traces!(
        cache::_RingCoronaCache,
        m::AbstractMetric,
        d::AbstractAccretionGeometry,
        β_angle;
    )

Calculates all of the base geodesics and optimizer geodesics for a given
point-approximation of a ring corona.
"""
function _ring_arm_traces!(
    cache::_RingCoronaCache,
    m::AbstractMetric{T},
    d::AbstractAccretionGeometry,
    β_angle;
    solver_opts...,
) where {T}
    _velfunc = rotatorfunctor(m, cache.x, cache.v, β_angle)

    # init a reusable integrator
    integ = _init_integrator(
        m,
        cache.x,
        _velfunc(0.0),
        d,
        # TODO: make this a parameter
        1_000_000.0,
        save_on = false,
        callback = domain_upper_hemisphere(),
        chart = chart_for_metric(m, 1_000_000.0),
        integrator_verbose = false,
        solver_opts...,
    )

    function _tracer(θ)
        init_v = _velfunc(θ)
        _solve_reinit!(integ, vcat(cache.x, init_v))
    end

    for i = 1:cache.N
        cache.gps[i] = _tracer(cache.angles[i])
    end
    # update the index
    cache._index[] = cache.N

    # do the minima and maxima extrema optimiser
    _golden_bracket!(cache, _tracer, Val{:minima}())
    _golden_bracket!(cache, _tracer, Val{:maxima}())
end

function unpack_traces(cache::_RingCoronaCache)
    last_index::Int = cache._index[]
    @views cache.angles[1:last_index], cache.gps[1:last_index]
end

function canonical_orders!(slices, rs)
    for I in slices
        sort!(I, by = i -> rs[i])
    end
    slices
end

# TODO: these are all incredibly allocating, when we could pre-allocate a
# buffer and reuse it
function _split_branches_further(inds, rs; cutoff_r = 1e3, delta = 0.01, min_length = 3)
    N = lastindex(inds)
    prev, first_i::Int = @views findmin(rs[inds])
    splits = Int[1, first_i]
    increasing::Bool = true

    for j = 0:N-1
        k = (j + first_i) % N + 1
        r = rs[inds[k]]

        if r > cutoff_r
            increasing = false
            prev = r
            continue
        end
        if increasing && (r < prev - delta)
            increasing = false
            push!(splits, k > 1 ? k - 1 : lastindex(inds))
        elseif !increasing && (r > prev + delta)
            increasing = true
            push!(splits, k > 1 ? k - 1 : lastindex(inds))
        end
        prev = r
    end

    sort!(splits)

    if length(splits) > 2
        slices = Vector{Int}[]
        i1 = first(splits)
        for i2 in splits[2:end]
            if i2 - i1 >= min_length
                if (N - i2) <= min_length
                    i2 = N
                end
                push!(slices, inds[i1:i2])
                i1 = i2
            end
        end
        if i1 < N
            push!(slices, inds[i1:end])
        end
        slices
    else
        [inds]
    end
end

function _split_arms_indices(angles, ρs)
    # split the sky into a left and right side, such that each side has strictly
    # sortable and monotonic r as a function of θ on the local sky

    # some values may be NaN so need to avoid those
    _start_index = findfirst(!isnan, ρs)

    _min_ρ_index = _start_index
    _max_ρ_index = _start_index
    _r_min = ρs[_start_index]
    _r_max = ρs[_start_index]

    for (i, r) in enumerate(ρs)
        r == NaN && continue
        if r < _r_min
            _min_ρ_index = i
            _r_min = r
        end
        if r > _r_max
            _max_ρ_index = i
            _r_max = r
        end
    end

    min_ρ_index, max_ρ_index =
        min(_min_ρ_index, _max_ρ_index), max(_min_ρ_index, _max_ρ_index)

    r1 = vcat(collect(max_ρ_index+1:lastindex(ρs)), collect(1:min_ρ_index-1))
    r2 = collect(min_ρ_index:max_ρ_index)

    l1 = canonical_orders!(_split_branches_further(r1, ρs), ρs)
    l2 = canonical_orders!(_split_branches_further(r2, ρs), ρs)
    vcat(l1, l2)
end

function _ring_arm!(
    cache::_RingCoronaCache,
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    β_angle;
    solver_opts...,
)
    _ring_arm_traces!(cache, m, d, β_angle; solver_opts...)

    all_angles, all_gps = unpack_traces(cache)

    # TODO: make the `_ring_arm_traces!` do this filtering automatically
    mask = [i.status == StatusCodes.IntersectedWithGeometry for i in all_gps]
    gps = all_gps[mask]
    angles = mod2pi.(all_angles[mask])

    I = sortperm(angles)
    angles = angles[I]
    gps = gps[I]

    ρs = map(i -> Gradus._equatorial_project(i.x), gps)

    # split the sky into a left and right side, such that each side has strictly
    # sortable and monotonic r as a function of θ on the local sky
    _, min_ρ_index = findmin(ρs)
    _, max_ρ_index = findmax(ρs)

    r1, r2 = _split_arms_indices(angles, ρs)

    # TODO: do we need views here? should go through this and work out where the memory gets copied
    left = @views Gradus._process_ring_traces(
        cache.setup,
        m,
        d,
        cache.v,
        gps[r1],
        ρs[r1],
        angles[r1],
    )
    right = @views Gradus._process_ring_traces(
        cache.setup,
        m,
        d,
        cache.v,
        gps[r2],
        ρs[r2],
        angles[r2],
    )

    Gradus.LongitudalArms(β_angle, left.r, left.t, left.ε, right.r, right.t, right.ε)
end

"""
    function corona_arms(
        setup::EmissivityProfileSetup{true},
        m::AbstractMetric,
        d::AbstractAccretionDisc,
        model::RingCorona,
        βs,
    )

Calculate [`LongitudalArms`](@ref) for a [`RingCorona`](@ref) for a given set
or range of ``\\beta`` angles. Here, ``\\beta`` is the angle relative to the
radial coordinate vector of the slice of geodesics being calculated (the
'slices of the orange' or 'beachball').

This function parallelises over CPU threads.
"""
function corona_arms(
    setup::EmissivityProfileSetup{true},
    m::AbstractMetric,
    d::AbstractAccretionDisc,
    model::RingCorona,
    βs;
    verbose = false,
    kwargs...,
)
    cache = _RingCoronaCache(setup, m, model)

    # copy one cache for each thread
    caches = [copy(cache) for i = 1:Threads.nthreads()-1]
    push!(caches, cache)

    progress_bar = init_progress_bar("β slices: ", length(βs), verbose)

    function _func(β)
        thread_cache = caches[Threads.threadid()]
        arm = _ring_arm!(thread_cache, m, d, β; kwargs...)
        thread_cache._index[] = 0
        ProgressMeter.next!(progress_bar)
        arm
    end

    arms = _threaded_map(_func, βs)
end

## BREAK

"""
    ang_diff(a1, a2)

Computes the angular difference between two angles using modulo arithmetic.
"""
function ang_diff(a1, a2)
    b1 = a1 > π / 2 ? π - (a1 % π) : a1
    b2 = a2 > π / 2 ? π - (a2 % π) : a2
    abs(b1 - b2)
end
"""
    PointSlice

A struct containing the raw values from a given ``\\beta`` slice trace.
"""
struct PointSlice{T}
    "Latitudal angles"
    θ::Vector{T}
    "Redshift"
    g::Vector{T}
    "Lorentz factor"
    γ::Vector{T}
    "Radial coordinate"
    r::Vector{T}
    # TODO: use AD to calculate δr instead of finite difference
    "Corona to disc time"
    t::Vector{T}
end

struct TimeDependentEmissivityBranch{T}
    "Latitudal angles"
    θ::Vector{T}
    "Radial coordinate"
    r::Vector{T}
    "Corona to disc time."
    t::Vector{T}
    "Emissivity"
    ε::Vector{T}
end

function _branch_emissivity(m::AbstractMetric, slice::PointSlice, I::Vector{Int}, spectrum)
    # TODO: this shares a lot of overlap code with other emissivity
    # implementations, and should be refactored to reuse more existing code
    map(eachindex(I)) do j
        j1, j2, j3, j4 = if j == 1
            1, 2, 1, 2
        elseif j != lastindex(I)
            j, j + 1, j, j - 1
        else
            j, j - 1, j, j - 1
        end

        # index conversion
        i = I[j]
        i1 = I[min(lastindex(I), j1)]
        i2 = I[min(lastindex(I), j2)]
        i3 = I[min(lastindex(I), j3)]
        i4 = I[min(lastindex(I), j4)]

        Δr = (abs(slice.r[i1] - slice.r[i2]) + abs(slice.r[i3] - slice.r[i4])) / 2
        weight =
            (ang_diff(slice.θ[i1], slice.θ[i2]) + ang_diff(slice.θ[i3], slice.θ[i4])) / 4

        A = _proper_area(m, slice.r[i] + Δr / 2, π / 2) * Δr
        weight * point_source_equatorial_disc_emissivity(
            spectrum,
            slice.θ[i],
            slice.g[i],
            A,
            slice.γ[i],
        )
    end
end

"""
    split_into_branches(m::AbstractMetric, slice::PointSlice, spectrum)

Split a given [`PointSlice`](@ref) into a number of
[`TimeDependentEmissivityBranch`](@ref) by cutting the curves of `(θ, r)` into
bijective branches.

Returns a vector of [`TimeDependentEmissivityBranch`](@ref).
"""
function split_into_branches(m::AbstractMetric, slice::PointSlice, spectrum)
    splits = _split_arms_indices(slice.θ, slice.r)
    map(splits) do I
        em = _branch_emissivity(m, slice, I, spectrum)
        TimeDependentEmissivityBranch(slice.θ[I], slice.r[I], slice.t[I], em)
    end
end

function empty_point_slice(; n = 1000, T = Float64)
    PointSlice(
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
    )
end

"""
    RingPointApproximationSlices

Used to hold all of the different [`PointSlice`](@ref) slices for their
respective values of ``\\beta``.
"""
struct RingPointApproximationSlices{T}
    "The slice angle."
    βs::Vector{T}
    "The slice of traces corresponding to a given β"
    traces::Vector{PointSlice{T}}
end

"""
    function interpolate_slice(
        pas::RingPointApproximationSlices{T},
        β;
        kwargs...
    )

Same as [`interpolate_slice!`](@ref) except will allocate the output for the
caller. The `kwargs` are forwarded to [`empty_point_slice`](@ref).
"""
function interpolate_slice(pas::RingPointApproximationSlices, β::Number; kwargs...)
    slice = empty_point_slice(; kwargs..., T = typeof(β))
    interpolate_slice!(slice, pas, β)
end

"""
    function interpolate_slice!(
        slice::PointSlice{T},
        pas::RingPointApproximationSlices{T},
        β,
    )

Interpolate a new slice into `slice` for ``\\beta`` given a set of calculated
slices in [`RingPointApproximationSlices`](@ref).
"""
function interpolate_slice!(
    slice::PointSlice{T},
    pas::RingPointApproximationSlices,
    β::T,
) where {T}
    idx = clamp(searchsortedlast(pas.βs, β), 1, length(pas.βs) - 1)
    x1, x2 = pas.βs[idx], pas.βs[idx+1]
    # interpolation weight
    w = (β - x1) / (x2 - x1)

    s1 = pas.traces[idx]
    s2 = pas.traces[idx+1]

    l1, u1 = extrema(s1.θ)
    l2, u2 = extrema(s2.θ)


    # interpolate the slice over the common support
    g1 = Gradus.NaNLinearInterpolator(s1.θ, s1.g, NaN)
    g2 = Gradus.NaNLinearInterpolator(s2.θ, s2.g, NaN)

    γ1 = Gradus.NaNLinearInterpolator(s1.θ, s1.γ, NaN)
    γ2 = Gradus.NaNLinearInterpolator(s2.θ, s2.γ, NaN)

    r1 = Gradus.NaNLinearInterpolator(s1.θ, s1.r, NaN)
    r2 = Gradus.NaNLinearInterpolator(s2.θ, s2.r, NaN)

    t1 = Gradus.NaNLinearInterpolator(s1.θ, s1.t, NaN)
    t2 = Gradus.NaNLinearInterpolator(s2.θ, s2.t, NaN)

    i::Int = 1
    for th in range(min(l1, l2), max(u1, u2), length(slice.θ))
        slice.θ[i] = th
        slice.g[i] = _linear_interpolate(g1(th), g2(th), w)
        slice.γ[i] = _linear_interpolate(γ1(th), γ2(th), w)
        slice.r[i] = _linear_interpolate(r1(th), r2(th), w)
        slice.t[i] = _linear_interpolate(t1(th), t2(th), w)
        i += 1
    end

    slice
end

function arrange_slice!(
    cache::_RingCoronaCache,
    m::AbstractMetric,
    d::AbstractAccretionDisc,
)
    _all_angles, _all_gps = unpack_traces(cache)

    # calculate the projected radial coordinate into the buffer
    for (i, gp) in enumerate(_all_gps)
        cache._buffer[i] = _equatorial_project(gp.x)
    end
    _all_rs = cache._buffer[1:length(_all_gps)]

    J = sortperm(_all_angles)
    unique!(i -> _all_rs[i], J)

    all_angles = _all_angles[J]
    all_rs = _all_rs[J]
    all_gps = _all_gps[J]

    # TODO: do we really need uniqueness here?

    # function for obtaining keplerian velocities
    _disc_velocity = _keplerian_velocity_projector(m, d)

    gammas = similar(all_rs)
    gs = map(enumerate(all_gps)) do dat
        i, p = dat
        v_disc = _disc_velocity(p.x)
        gammas[i] = lorentz_factor(m, p.x, v_disc)
        energy_ratio(m, p, cache.v, v_disc)
    end

    mask = [i.status == StatusCodes.IntersectedWithGeometry for i in all_gps]
    ts = [i.x[1] for i in all_gps]

    # mask the traces that didn't hit the disc
    # but not on angles, since we use those for interpolating
    @. gs[!mask] = all_rs[!mask] = ts[!mask] = gammas[!mask] = NaN
    PointSlice(all_angles, gs, gammas, all_rs, ts)
end

"""
    function corona_slices(
        setup::EmissivityProfileSetup{true},
        m::AbstractMetric,
        d::AbstractAccretionDisc,
        model::RingCorona,
        βs,
    )

Calculates a [`RingPointApproximationSlices`](@ref) for a [`RingCorona`](@ref) for a
given set or range of ``\\beta`` angles. Here, ``\\beta`` is the angle relative
to the radial coordinate vector of the slice of geodesics being calculated (the
'slices of the orange' or 'beachball').

This function parallelises over CPU threads.
"""
function corona_slices(
    setup::EmissivityProfileSetup{true},
    m::AbstractMetric,
    d::AbstractAccretionDisc,
    model::RingCorona,
    βs;
    verbose = false,
    kwargs...,
)
    cache = _RingCoronaCache(setup, m, model)
    # copy one cache for each thread
    caches = [copy(cache) for i = 1:Threads.nthreads()-1]
    push!(caches, cache)

    progress_bar = init_progress_bar("β slices: ", length(βs), verbose)
    function _func(β)
        thread_cache = caches[Threads.threadid()]
        _ring_arm_traces!(thread_cache, m, d, β; kwargs...)

        slice = @views arrange_slice!(thread_cache, m, d)

        # reset for next iteration
        thread_cache._index[] = 0
        ProgressMeter.next!(progress_bar)

        slice
    end

    traces = Gradus._threaded_map(_func, βs)
    RingPointApproximationSlices(βs, traces)
end


struct RingApproximation{T}
    βs::Vector{T}
    branches::Vector{Vector{TimeDependentEmissivityBranch{T}}}
end

function emissivity_at(ra::RingApproximation, r)
    sum(enumerate(ra.branches)) do dat
        i, branches = dat
        i1, i2 = if i >= (length(ra.branches) - 1)
            (i - 1, i)
        else
            (i, i + 1)
        end

        v = sum(branches) do branch
            first_i = findfirst(!isnan, branch.r)
            last_i = findlast(!isnan, branch.r)
            if !isnothing(first_i) &&
               !isnothing(last_i) &&
               (r >= branch.r[first_i]) &&
               (r <= branch.r[last_i])
                Gradus._make_interpolation(branch.r, branch.ε)(r)
            else
                zero(typeof(r))
            end
        end

        v
    end
end

function make_approximation(
    m::AbstractMetric,
    slices::RingPointApproximationSlices,
    spectrum;
    βs = collect(range(extrema(slices.βs)..., 1000)),
)
    slice = empty_point_slice()
    branches = map(βs) do beta
        interpolate_slice!(slice, slices, beta)
        split_into_branches(m, slice, spectrum)
    end
    RingApproximation(βs, branches)
end
