function _find_offset_for_measure(
    measure,
    m::AbstractMetric,
    x,
    d::AbstractAccretionGeometry,
    θₒ;
    zero_atol = 1e-7,
    offset_max = 20.0,
    initial_r = offset_max / 2,
    max_time = 2 * x[2],
    β₀ = 0,
    α₀ = 0,
    μ = 0,
    solver_opts...,
)
    function _velfunc(r::T) where {T}
        α, β = _rθ_to_αβ(r, θₒ; α₀ = α₀, β₀ = β₀)
        # need constrain_all here since `_solve_reinit!` doesn't normalize
        constrain_all(m, x, map_impact_parameters(m, x, α, β), μ)
    end

    # init a reusable integrator
    integ =
        _init_integrator(m, x, _velfunc(0.0), d, max_time; save_on = false, solver_opts...)

    function f(r)
        if r < 0
            return -1000 * r
        end
        gp = _solve_reinit!(integ, vcat(x, _velfunc(r)))
        measure(gp)
    end

    # use adaptive Order0 method : https://juliamath.github.io/Roots.jl/dev/reference/#Roots.Order0
    r0 = Roots.find_zero(f, initial_r, Roots.Order0(); atol = zero_atol)

    gp0 = _solve_reinit!(integ, vcat(x, _velfunc(r0)))
    r0, gp0
end

function _find_offset_for_radius(
    m::AbstractMetric,
    x,
    d::AbstractAccretionGeometry,
    rₑ,
    θₒ;
    warn = true,
    kwargs...,
)
    function _measure(gp::GeodesicPoint{T}) where {T}
        r = if gp.status == StatusCodes.IntersectedWithGeometry
            _equatorial_project(gp.x)
        else
            zero(T)
        end
        rₑ - r
    end
    r0, gp0 = _find_offset_for_measure(_measure, m, x, d, θₒ; kwargs...)

    if warn && (r0 < 0)
        @warn("Root finder found negative radius for rₑ = $rₑ, θₑ = $θₒ")
    end
    if !isapprox(_measure(gp0), 0.0, atol = 1e-5)
        warn && @warn("Poor offset radius found for rₑ = $rₑ, θₑ = $θₒ")
        return NaN, gp0
    end
    r0, gp0
end

"""
    find_offset_for_radius(
        m::AbstractMetric,
        x,
        d::AbstractAccretionDisc,
        radius,
        θₒ;
        zero_atol = 1e-7,
        offset_max = 20.0,
        kwargs...,
    )

Find the offset ``r_\\text{o}`` on the observer's image plane that gives impact parameters
```math
\\alpha = r_\\text{o} \\cos \\theta_\\text{o},
\\quad \\text{and} \\quad
\\beta = r_\\text{o} \\sin \\theta_\\text{o},
```
that trace a geodesic to intersect the disc geometry at an emission radius ``r_\\text{e}``.

Returns ``NaN`` if no offset exists. This may occur of the geometry is not present at this location
(though this more commonly gives a bracketing interval error), or if the emission radius is within 
the event horizon.
"""
find_offset_for_radius(
    m::AbstractMetric,
    x,
    d::AbstractAccretionGeometry,
    rₑ,
    θₒ;
    kwargs...,
) = _find_offset_for_radius(m, x, d, rₑ, θₒ; kwargs...)
find_offset_for_radius(
    m::AbstractMetric,
    x,
    d::AbstractThickAccretionDisc,
    rₑ,
    θₒ;
    kwargs...,
) = _find_offset_for_radius(m, x, datumplane(d, rₑ), rₑ, θₒ; kwargs...)

"""
    impact_parameters_for_radius!(
        αs::AbstractVector,
        βs::AbstractVector,
        m::AbstractMetric,
        x::AbstractVector{T},
        d::AbstractAccretionDisc,
        rₑ;
        kwargs...,
    )

Finds ``\\alpha`` and ``\\beta`` impact parameters which trace geodesics to a given emission radius
`rₑ` on the accretion disc.

This function assigns the values in-place to `αs` and `βs`.
The keyword arguments are forwarded to [`find_offset_for_radius`](@ref).

This function is threaded.
"""
function impact_parameters_for_radius!(
    αs::AbstractVector,
    βs::AbstractVector,
    m::AbstractMetric,
    x::AbstractVector{T},
    d::AbstractAccretionDisc,
    rₑ;
    β₀ = zero(T),
    α₀ = zero(T),
    kwargs...,
) where {T}
    if size(αs) != size(βs)
        throw(DimensionMismatch("α, β must have the same dimensions and size."))
    end
    θs = range(0, 2π, length(αs))
    @inbounds @threads for i in eachindex(θs)
        θ = θs[i]
        r, _ = find_offset_for_radius(m, x, d, rₑ, θ; α₀ = α₀, β₀ = β₀, kwargs...)
        αs[i] = r * cos(θ) + α₀
        βs[i] = r * sin(θ) + β₀
    end
    (αs, βs)
end

"""
    function impact_parameters_for_radius(
        m::AbstractMetric,
        x::AbstractVector{T},
        d::AbstractAccretionDisc,
        radius;
        N = 500,
        kwargs...,
    )

Pre-allocating version of [`impact_parameters_for_radius!`](@ref), which allocates
for `N` impact paramter pairs. Filters for `NaN` and returns `(αs, βs)`.
"""
function impact_parameters_for_radius(
    m::AbstractMetric,
    x::AbstractVector{T},
    d::AbstractAccretionDisc,
    radius;
    N = 500,
    kwargs...,
) where {T}
    α = zeros(T, N)
    β = zeros(T, N)
    impact_parameters_for_radius!(α, β, m, x, d, radius; kwargs...)
    α, β
end

function impact_parameters_for_radius_obscured(
    m::AbstractMetric,
    x::AbstractVector{T},
    d::AbstractAccretionDisc,
    radius;
    N = 500,
    β₀ = zero(T),
    α₀ = zero(T),
    kwargs...,
) where {T}
    α = zeros(T, N)
    β = zeros(T, N)
    impact_parameters_for_radius!(α, β, m, x, d, radius; β₀ = β₀, α₀ = α₀, kwargs...)

    gps = tracegeodesics(
        m,
        x,
        i -> map_impact_parameters(m, x, α[i], β[i]),
        datumplane(d, radius),
        2x[2];
        ensemble = EnsembleEndpointThreads(),
        save_on = false,
        trajectories = length(α),
    )
    n = _cartesian_surface_normal(d, radius)
    I = [!_is_visible(m, d, gp, n) for gp in gps]
    α[I] .= NaN
    β[I] .= NaN
    α, β
end

"""
If dot product between surface normal and the final velocity vector is positive, 
the point is visible.
"""
function _is_visible(m::AbstractMetric, d, gp::AbstractGeodesicPoint, n::SVector{3})
    gp_new =
        tracegeodesics(m, gp.x_init, gp.v_init, d, gp.λ_max; save_on = false) |>
        unpack_solution
    # geodesics sufficiently close, test the angle
    # v = SVector(gp_new.v[2], gp_new.v[3], gp_new.v[4]) ./ gp_new.v[1]
    # v_geodesic = _spher_to_cart_jacobian(gp_new.x[3], gp_new.x[4], gp_new.x[2]) * v
    # n_rot = _rotate_about_spinaxis(n, gp_new.x[4])

    X1 = to_cartesian(gp.x)
    X2 = to_cartesian(gp_new.x)
    dist = sum(i -> i^2, X1 .- X2)

    # return _fast_dot(v_geodesic, n_rot) < 0 && 
    return isapprox(dist, 0, atol = 1e-12)
end

function jacobian_∂αβ_∂gr(
    m::AbstractMetric,
    x::SVector{4},
    d::AbstractAccretionDisc,
    α,
    β,
    max_time;
    redshift_pf = ConstPointFunctions.redshift(m, x),
    solver_opts...,
)
    domain_limiter = d isa AbstractThickAccretionDisc ? domain_upper_hemisphere() : nothing

    # these type hints are crucial for forward diff to be type stable
    function _jacobian_f(impact_params::SVector{2,T})::SVector{2,T} where {T}
        sol = tracegeodesics(
            m,
            x,
            map_impact_parameters(m, x, impact_params[1], impact_params[2]),
            d,
            max_time;
            save_on = false,
            callback = domain_limiter,
            solver_opts...,
        )
        gp = unpack_solution(sol)
        ρ = _equatorial_project(gp.x)

        # if within the inner radius of the disc, the redshift might be more or 
        # less undefined, or at the very least discontinuous, causing havoc with
        # the gradients. furthermore, we shouldn't expect to have any geodesics that
        # fall within the inner radius thanks to the root finder. so what we do
        # instead is set δg to zero, which causes the inverse jacobian to be infinite
        # which we later use as a check for "visibility"
        # 
        # todo: this is not perfect, but works well enough
        g = if (d isa AbstractThickAccretionDisc) && (ρ < inner_radius(d))
            zero(T)
        else
            redshift_pf(m, gp, max_time)
        end

        SVector(ρ, g)
    end

    # compute |∂(r,g) / ∂(α,β)|⁻¹
    j_matrix = ForwardDiff.jacobian(_jacobian_f, SVector(α, β))
    abs(inv(det(j_matrix)))
end

function _make_target_objective(
    target::SVector,
    m::AbstractMetric,
    x0,
    args...;
    d_tol = 1e-2,
    max_time = 2x0[2],
    callback = nothing,
    μ = 0.0,
    solver_opts...,
)
    # convenience velocity function
    velfunc = (α, β) -> begin
        constrain_all(m, x0, map_impact_parameters(m, x0, α, β), μ)
    end

    # used to track how close the current solver got
    closest_approach = Ref(x0[2])
    # convert target to cartesian once
    target_cart = to_cartesian(target)
    distance_callback = ContinuousCallback(
        (u, λ, integrator) -> begin
            # get vector distance between current x and target
            k = target_cart - to_cartesian(u)
            distance = √sum(i -> i^2, k)
            # update best
            closest_approach[] = min(closest_approach[], distance)
            # if within d_tol, just immediately terminate
            distance - d_tol
        end,
        terminate!,
        interp_points = 8,
        save_positions = (false, false),
    )

    # init a reusable integrator
    integ = _init_integrator(
        m,
        x0,
        velfunc(0.0, 0.0),
        args...,
        (0.0, max_time);
        save_on = false,
        callback = merge_callbacks(callback, distance_callback),
        solver_opts...,
    )

    function _solver(α, β)
        v = velfunc(α, β)
        _solve_reinit!(integ, vcat(x0, v))
    end

    # map impact parameters to a closest approach 
    function _target_objective(impact_params)
        α = impact_params[1]
        β = impact_params[2]
        # reset the closest approach
        closest_approach[] = x0[2]
        _ = _solver(α, β)
        closest_approach[]
    end

    return _target_objective, _solver
end

function optimize_for_target(
    target::SVector,
    m::AbstractMetric,
    x0::SVector{4,T},
    args...;
    optimizer = NelderMead(),
    p0 = zeros(T, 2),
    kwargs...,
) where {T}
    f, solver = _make_target_objective(target, m, x0, args...; kwargs...)
    res = optimize(f, p0, optimizer)
    out = Optim.minimizer(res)
    # return α, β, accuracy
    α, β = out[1], out[2]
    accuracy = Optim.minimum(res)

    α, β, unpack_solution(solver(α, β)), accuracy
end

function impact_parameters_for_target(
    target::SVector,
    m::AbstractMetric,
    x0,
    args...;
    kwargs...,
)
    α, β, _, accuracy = optimize_for_target(target, m, x0, args...; kwargs...)
    α, β, accuracy
end

export find_offset_for_radius, impact_parameters_for_radius, impact_parameters_for_radius!
