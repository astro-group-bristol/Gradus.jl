const EXTREMIZER_ITERS = 90
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
    "Used to track which angles and gps we have written so far. `_index == 0` means the array is empty."
    _index::Ref{Int}
    "Number of geodesics for all the arms"
    N::Int
end

function Base.copy(cache::_RingCoronaCache)
    _RingCoronaCache(
        cache.setup,
        cache.x,
        cache.v,
        deepcopy(cache.angles),
        deepcopy(cache.gps),
        Ref(cache._index[]),
        cache.N,
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
)
    x, v = sample_position_velocity(m, model)

    # TODO: infer this
    GPType = GeodesicPoint{eltype(x),Nothing}

    N = (setup.n_samples - 2 * EXTREMIZER_ITERS)
    if (N <= 0)
        error(
            "Too few samples for selected method. Pass `n_samples` of more than $(2 * EXTREMIZER_ITERS)",
        )
    end
    if !iseven(N)
        error("`n_samples` must be an even number for this method (given $N)")
    end

    angles = zeros(eltype(x), setup.n_samples)
    gps = Vector{GPType}(undef, setup.n_samples)

    num_per_arm = div(N, 2)

    # setup the angles we are going to sample evenly over the sky
    left_arm_angles = range(h, π - h, num_per_arm)
    for (i, a) in enumerate(left_arm_angles)
        angles[i] = a
        angles[i+num_per_arm] = mod2pi(a + π)
    end

    _RingCoronaCache(setup, x, v, angles, gps, Ref(0), N)
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

    # TODO: this should be a seperate function
    # find the true index from the sample index
    sample_gps = (
        i for i in @views(cache.gps[1:cache.N]) if
        i.status == StatusCodes.IntersectedWithGeometry
    )
    best_estimate = _equatorial_project(first(sample_gps).x)
    best_angle = first(cache.angles)
    for i = 1:cache.N
        gp = cache.gps[i]
        if gp.status == StatusCodes.IntersectedWithGeometry
            rho = _equatorial_project(gp.x)
            if comparator(rho, best_estimate)
                best_estimate = rho
                best_angle = cache.angles[i]
            end
        end
    end

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

    delta = π / min(cache.N, 60)

    # take a bracketing interval
    a = best_angle - delta
    b = best_angle + delta
    c_value = 0.0

    n = 0
    iters = 0
    too_many_iters = false
    while (n < EXTREMIZER_ITERS)
        c = b - (b - a) * INVPHI
        d = a + (b - a) * INVPHI

        c_gp, c_value = _objective(c)
        d_gp, d_value = _objective(d)

        if comparator(c_value, d_value)
            if (c_value != best_estimate)
                n += 1
                _add_to_cache!(cache, c, c_gp)
            elseif too_many_iters
                n += 1
            end
            b = d
        else
            if (d_value != best_estimate)
                n += 1
                _add_to_cache!(cache, d, d_gp)
            elseif too_many_iters
                n += 1
            end
            a = c
        end

        iters += 1
        if !too_many_iters && (iters > EXTREMIZER_ITERS)
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

    min_ρ_index, max_ρ_index = min(min_ρ_index, max_ρ_index), max(min_ρ_index, max_ρ_index)

    # TODO: do we need views here? should go through this and work out where the memory gets copied
    left = @views Gradus._process_ring_traces(
        cache.setup,
        m,
        d,
        cache.v,
        vcat(gps[1:(min_ρ_index-1)], gps[max_ρ_index+1:end]),
        vcat(ρs[1:(min_ρ_index-1)], ρs[max_ρ_index+1:end]),
        vcat(angles[1:(min_ρ_index-1)], angles[max_ρ_index+1:end]),
    )
    right = @views Gradus._process_ring_traces(
        cache.setup,
        m,
        d,
        cache.v,
        gps[min_ρ_index:max_ρ_index],
        ρs[min_ρ_index:max_ρ_index],
        angles[min_ρ_index:max_ρ_index],
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

    function _func(β)
        thread_cache = caches[Threads.threadid()]
        arm = _ring_arm!(thread_cache, m, d, β)
        thread_cache._index[] = 0
        arm
    end

    arms = _threaded_map(_func, βs)
end
