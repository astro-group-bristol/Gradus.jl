"""
    find_offset_for_radius(
        m::AbstractMetricParams,
        u,
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
function find_offset_for_radius(
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionDisc,
    rₑ,
    θₒ;
    zero_atol = 1e-7,
    offset_max = 20.0,
    max_time = 2 * u[2],
    μ = 0.0,
    solver_opts...,
)
    measure = gp -> begin
        r = if gp.status == StatusCodes.IntersectedWithGeometry
            gp.u2[2]
        else
            0.0
        end
        rₑ - r
    end

    velfunc = r -> begin
        α = r * cos(θₒ)
        β = r * sin(θₒ)
        constrain_all(m, u, map_impact_parameters(m, u, α, β), μ)
    end
    # init a reusable integrator
    integ = _init_integrator(
        m,
        u,
        velfunc(0.0),
        d,
        (0.0, max_time);
        save_on = false,
        solver_opts...,
    )

    f = r -> begin
        gp = _solve_reinit!(integ, vcat(u, velfunc(r)))
        measure(gp)
    end

    # use adaptive Order0 method : https://juliamath.github.io/Roots.jl/dev/reference/#Roots.Order0
    r0 = Roots.find_zero(f, offset_max / 2; atol = zero_atol)

    gp0 = _solve_reinit!(integ, vcat(u, velfunc(r0)))
    if !isapprox(measure(gp0), 0.0, atol = 10 * zero_atol)
        return NaN, gp0
    end
    r0, gp0
end

"""
    impact_parameters_for_radius!(
        αs::AbstractVector,
        βs::AbstractVector,
        m::AbstractMetricParams,
        u::AbstractVector{T},
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
    m::AbstractMetricParams,
    u::AbstractVector{T},
    d::AbstractAccretionDisc,
    rₑ;
    kwargs...,
) where {T}
    if size(αs) != size(βs)
        throw(DimensionMismatch("α, β must have the same dimensions and size."))
    end
    θs = range(0, 2π, length(αs))
    @inbounds @threads for i in eachindex(θs)
        θ = θs[i]
        r, _ = find_offset_for_radius(m, u, d, rₑ, θ; kwargs...)
        αs[i] = r * cos(θ)
        βs[i] = r * sin(θ)
    end
    (αs, βs)
end

"""
    function impact_parameters_for_radius(
        m::AbstractMetricParams,
        u::AbstractVector{T},
        d::AbstractAccretionDisc,
        radius;
        N = 500,
        kwargs...,
    )

Pre-allocating version of [`impact_parameters_for_radius!`](@ref), which allocates
for `N` impact paramter pairs. Filters for `NaN` and returns `(αs, βs)`.
"""
function impact_parameters_for_radius(
    m::AbstractMetricParams,
    u::AbstractVector{T},
    d::AbstractAccretionDisc,
    radius;
    N = 500,
    kwargs...,
) where {T}
    α = zeros(T, N)
    β = zeros(T, N)
    impact_parameters_for_radius!(α, β, m, u, d, radius; kwargs...)
    (filter(!isnan, α), filter(!isnan, β))
end

function jacobian_∂αβ_∂gr(
    m,
    u,
    d,
    α,
    β,
    max_time;
    diff_order = 5,
    μ = 0.0,
    redshift_pf = ConstPointFunctions.redshift(m, u),
    solver_opts...,
)
    velfunc = (α, β) -> begin
        constrain_all(m, u, map_impact_parameters(m, u, α, β), μ)
    end

    # init a reusable integrator
    integ = _init_integrator(
        m,
        u,
        velfunc(0.0, 0.0),
        d,
        (0.0, max_time);
        save_on = false,
        solver_opts...,
    )

    # map impact parameters to r, g
    f = ((α, β),) -> begin
        v = velfunc(α, β)
        gp = _solve_reinit!(integ, vcat(u, v))
        g = redshift_pf(m, gp, max_time)
        # return r and g*
        @SVector [gp.u2[2], g]
    end

    # choice between FiniteDifferences and FiniteDiff is tricky
    # since FiniteDiff is so much faster, but seems to yield really bad jacobians
    # for this specific problem, so instead stenciling with a given order
    cfdm = FiniteDifferences.central_fdm(diff_order, 1)
    j = FiniteDifferences.jacobian(cfdm, f, @SVector([α, β])) |> first
    abs(inv(det(j)))
end

function make_target_objective(
    target::SVector,
    m::AbstractMetricParams,
    u0,
    args...;
    d_tol = 1e-2,
    max_time = 2u0[2],
    callback = nothing,
    μ = 0.0,
    solver_opts...,
)
    # convenience velocity function
    velfunc = (α, β) -> begin
        constrain_all(m, u0, map_impact_parameters(m, u0, α, β), μ)
    end

    # used to track how close the current solver got
    closest_approach = Ref(u0[2])
    # convert target to cartesian once
    target_cart = to_cartesian(target)
    distance_callback = ContinuousCallback(
        (u, λ, integrator) -> begin
            # get vector distance between current u and target
            k = target_cart - to_cartesian(u)
            distance = √sum(i -> i^2, k)
            # update best
            closest_approach[] = min(closest_approach[], distance)
            # if within d_tol, just immediately terminate
            distance - d_tol
        end,
        terminate!,
        interp_points = 8,
        save_positions = (true, false),
    )

    # init a reusable integrator
    integ = _init_integrator(
        m,
        u0,
        velfunc(0.0, 0.0),
        args...,
        (0.0, max_time);
        save_on = false,
        callback = merge_callbacks(callback, distance_callback),
        solver_opts...,
    )

    # map impact parameters to a closest approach 
    f = ((α, β),) -> begin
        # reset the closest approach
        closest_approach[] = u0[2]
        v = velfunc(α, β)
        _ = _solve_reinit!(integ, vcat(u0, v))
        closest_approach[]
    end

    return f
end

function impact_parameters_for_target(
    target::SVector,
    m::AbstractMetricParams,
    u0,
    args...;
    optimizer = NelderMead(),
    kwargs...,
)
    f = make_target_objective(target, m, u0, args...; kwargs...)
    res = optimize(f, [0.0, 0.0], optimizer)
    out = Optim.minimizer(res)
    # return α, β, accuracy
    out[1], out[2], Optim.minimum(res)
end

export find_offset_for_radius, impact_parameters_for_radius, impact_parameters_for_radius!
